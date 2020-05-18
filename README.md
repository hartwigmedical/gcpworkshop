## GCP Workshop 

### Using the console

Your first entry point to GCP will likely be the console. It gives you a tool to both interact and monitor various GCP services and even
create a quick shell for more advances operations.

Hit the following link to start up the console: https://console.cloud.google.com/

### Projects, Users, Roles and Service Accounts

When you begin in the console you'll notice that you are working within a *Project*. Projects group together services and users along with 
billing. Projects are a handy way to organize resources around groups of users with similar permissions and goals. We recommend taking a
modular approach when it comes to projects, many small cohesive ones over a single monolith. Good to note that a user can belong to many 
projects.

IAM policies are used to organize users permissions. It uses the concepts of roles, which group permissions into cohesive abilities to 
perform tasks. There are many roles, and you can even define custom ones when you don't have the right combination out of the box. An IAM 
policy can be applied to any resources: a project; a VM; a network; a bucket; etc.
[This page](https://cloud.google.com/iam/docs/understanding-roles) gives the overview of how roles map to permissions.

Service Accounts are an important concept and become useful when you want automate and scale. A service account can be used as a proxy to 
your end user credentials, to use GCP services via a private key rather than an explicit login. This has a variety of applications, from
daisy chaining services to perform more complex task (ie a VM which accesses storage) or to launch hundreds of VMs without the need for 
human interaction.

### Using the command line

The fastest way to get to a shell is to use the "cloud shell" functionality. The cloud shell already has all tools installed and will 
default to the project you have open. It works by provisioning a small VM.

![Cloud Shell](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/cloud-shell.png)

Open a cloud shell and try the following:

```
gcloud auth list
gsutil ls
```

That said, the shell still runs in the browser so can be a little slow and not suited to long running interactions. You can also do 
everything from a terminal on your own maching with the `gcloud` and `gsutil` tools.

Open a terminal and try the following:
```
gcloud auth login
gsutil ls
```

### Google Compute Engine (GCE)

GCE gives you the ability to create VMs, make images, and deploy containers. Next we'll create our own VM and learn how to access it via 
SSH.

Start with reviewing all the options from the Create VM on the [console](https://console.cloud.google.com/compute/instancesAdd):

![Create a VM](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/compute-demo-1.png)

Important Concepts:
* VMs also need to deployed to a region, but also a zone, which is part of the datacenter
* There are many pre-defined machine types, all pricing is linear to CPU and memory. You can define custom machines via the command line.
* You can define a disk when you define your VM, but this can also be done independently, and read-only disks can be shared.
* Inside the VM, service interactions are performed by the project service account. A service account can run headless, but otherwise has
permissions just like a normal user. 

Let's create one now with a 100GB disk.

You can also create and manage vms via the command line. 

```shell script
gcloud compute instances create {yourname}-instance --zone=europe-west4-a --boot-disk-size 100 --machine-type=n1-standard-2
```

The [VM instance overview](https://console.cloud.google.com/compute/instances) in the console is a nice way to see the status of all VMs:

![VM Instances](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/compute-demo-2.png)

The SSH button gives you a quick browser shell to the machine. Again, the browser has its limitations, so sometimes it's better to tunnel 
into your VM from your own terminal:

```shell script
gcloud compute --project "nki-atlas" ssh --zone "europe-west4-a" "{yourname}-vm"
```

Try SSH'ing into your new VM via the terminal with the gcloud command.

#### Virtual Private Cloud Networking

When working with sensitive data you most likely do not want your VM exposed anywhere on the public internet. ![VPC networking](https://cloud.google.com/vpc) 
provides you the ability to create a completely private network for all your VMs. We won't go into the details today, but the simplest way
you can protect your vms is to disable its public IP. 

```shell script
gcloud compute --project "nki-atlas" ssh --zone "europe-west4-a" --no-address "{yourname}-vm" 
```

Or via the Network tab in the console instance creation page.

#### Install Miniconda and Samtools 
We need to initialize the system so that we have access to bioinformatics tools, therefore we will install Miniconda on the VM. Walk through the instructions for installing and adding the updating you PATH information.
```shell script
# Copy Miniconda installer from the workshop bucket & install
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install Miniconda
## 1. Hit [Enter] to begin
## 2. Hit 'q' to skip to bottom and type 'yes' 
## 3. Hit [Enter] to confirm installation path
## 4. Type 'yes' to run initializations 
bash Miniconda3-latest-Linux-x86_64.sh

# After installation re-initialize your PATH
source .bashrc

# Set up conda environment
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install samtools using conda
conda install -y samtools=1.10
conda install -y fastqc
```

### Images

Images are a handy way to save and share state of a VM. You can take an image on the command line by:

```shell script
 gcloud compute images create {yourname}-image-1 --family={yourname}-image --source-disk={yourname}-vm --source-disk-zone=europe-west4 --storage-location=europe-west4
```

### GCE Cost Savings

During our migration to GCP, we found two major ways to cut back on our compute costs.

#### Pre-emptible VMS

The single biggest win we had was moving to all pre-emptible VMs. A pre-emptible VM, is one that GCP can claim back at any time, and will
always be shut-down at 24 hours. This is how Google sells excess capacity, and claims it back when they need it. It is a great fit for 
analysis workloads running between 1-24 hours, and the savings are big (80%). 

You can mark a VM pre-emptible on VM creation *Create VM -> Management -> Availability Policy* and through the `glcoud` CLI.

![Make a pre-emptible VM](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/compute-demo-3.png)

While creating them is easy, its worth having some automation around them to manage pre-emptions. In our pipline, we handle the pre-empted
signal by polling the GCE API, then re-starting the workload in a new zone.

#### Local SSDs

Local SSDs give you the best disk performance at the lowest cost. That said, you need to perform some steps if you need a disk larger than 
375GB and data will be completely transient. A VM with local SSDs cannot be restarted and all data there is lost. 

They are a great way to both speed up and reduce cost of a completely transient workload.

### Google Cloud Storage (GCS)

Storage in the cloud is one of its most powerful tools and concepts. GCS is an object storage system which works with buckets, paths and files.

#### Creating a Bucket
The first time, its best to create a bucket via the console:

![Create a bucket](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/storage-demo-1.png)

Important concepts:
* The region is where your bucket will be created. Best to pick a local region.
* The storage class impacts how your data is stored. Standard will be used for frequently accessed data. Nearline and Coldine for less 
accessed.
* The access control can be per bucket or per object, more on that later/
* By default data in encrypted with Google's internal keys. Customer Managed Encryption is also available, which we leverage to ensure no 
one can access our other than those we enable on the ACL. 

Create a standard storage bucket called `{yourname}-gcpdemo` in europe-west4, with object level access control and Google encryption. 
Bucket names have some restrictions: they must be globally unique; they should not have underscores; they must be lower-case. 
This is because you can expose data via direct URL access and must conform to the rules of DNS.

From there its best to move to the command line. After you've create your bucket, you can use the `gsutil` tool to upload, download, and
copy data between buckets.  

```
touch myfile.txt
gsutil -m cp myfile.txt gs://{yourname}-gcpdemo/
gsutil ls gs://{yourname}-gcpdemo/
```

The `-m` flag above will do the download in parallel. Keep in mind that this will be heavy on the bandwidth, but have a huge impact on the
speed of your operation. You can also get things to go even faster by uploading data as a composite object (aka multipart). The following
command will upload any file over 150M as a composite object, which also allows it to be downloaded as such. **Note**: You will not have an
MD5 checksum if a file is uploaded as composite.

```shell script
gsutil -m -o GSUtil:parallel_composite_upload_threshold=150M cp myfile.txt gs://{yourname}-gcpdemo/
```

You can also browse your data from within the console. Find Storage in your console and you can see your new bucket and file there.

Access to a bucket can controlled how we've done it or as an IAM policy on the bucket. Try the following commands to see the difference:

```shell script
# Object ACLs on your new bucket
gsutil acl get gs://{yourname}-gcpdemo/myfile.txt
# Bucket policy to our demo bucket 
gsutil -u nki-atlas iam get gs://nki-demo-data
```

Object ACLs have 3 permissions: READER, WRITER and OWNER. IAM policies can be any number of role bindings, giving much more flexibility and
specificity in access. 

### Accessing HMF data

HMF's raw data is stored in GCS in CRAM format. We give requesters permissions directly on the data rather than making copies to prevent
costly duplication. Currently we grant these permissions via ACLs as our buckets contain more than one patient's data. 

In this final exercise, we'll download some HMF cell line data to a VM, perform an operation, then copy our results back to GCS. We see this
as a common pattern, as it avoids any egress cost, is fast, and can be scaled horizontally.

Find the image you created earlier in the console and create a new VM instance based on it (hint, the console offers a quick way of doing
this from the instances page).  

To access HMF data for the demo you'll have to login as your own account in the VM. Do a quick list of authorized users:

```shell script
gcloud auth list
```

By default your VM will use a service account. You can apply to allow a service account access to your data request given that it is 
only accessible to users listed in the agreement. Run the login command to authenticate as yourself.

```shell script
gcloud auth login
```

HMF allows you to access the data with an approved service account as well. An approved service account should only be accessible by the
user 

When a data request is approved, we create a bucket for you with a file called `manifest.json`. The manifest contains a list of all files 
for which the requester has been granted access as GCS urls. Download the manifest for this workshop like so:

```shell script
gsutil -u nki-atlas cp gs://nki-demo-data/manifest.json ./
```

You'll notice in this `gsutil` command we used the `-u nki-atlas`. Data we provide externally are stored *requester pays* buckets. With a 
requester pays bucket the requester will be billed for any egress costs (0.12 GB), and this argument specifies the billing project. There
is no egress cost to download the data to a VM in the Netherlands, but you still need to specify the billing project. 

View the content of the JSON file.
```shell script
cat manifest.json
```

There you'll see the GCS urls mentioned earlier. Copy the url for the reference cram and 
the corresponding crai and copy them down to the VM with `gsutil`. See if you can create the command yourself (tip you'll definitely want
to use `-m` on this one!)

**Note:** If you do not need access to the raw CRAM files and need only the processed genomics data, you can download the `manifest.json` file directly and copy the url to download the necessary folder instead of placing the data first in a VM. 

In this exercise we want to slice out a region of the mapped reads, specifically the HLA regions. Therefore we can run the following samtools 
command to extract a particular set of read coordinates. 
```shell script
# Get mapping statistics
samtools index COLO829R.cram
samtools flagstat COLO829R.cram

# Next slice out the reads for an HLA region on chromosome 6.
samtools view -b COLO829R.cram 6:29909037-29913661 > COLO829R.sliced.bam
```

After extracting the alignments we can run a quality analysis of the output using fastqc.
```shell script
fastqc -o ./ -f bam COLO829R.sliced.bam
```

With access to the Anaconda package repository, there are many different analyses you can perform on the data. Next we need to upload 
the output data to our bucket so we can view and use it locally. 
```shell script
# Copy the sliced BAM data to your bucket
gsutil -m cp COLO829R.sliced.bam gs://{yourname}-gcpdemo

# Copy the FastQC report to your bucket
gsutil -m cp COLO829R.sliced.cram_fastqc.html gs://{yourname}-gcpdemo
```
That's it! Your data is now saved to GCS so you can stop or delete the VM. 

### Stop the VM

Navigate back to the console and go to Compute Engine. Select your VM and stop it from the toolbar. You'll noticed the VM is now stopped
but still in the list. *This VM still accrues storage cost*. The disk is still stored such that when you start it again, all the state is
maintained.

### Start a Jupyter Notebook instance
A more advanced exercise is to start a VM with Jupyter Notebook that can be modified for fit you computational requirements. Setup is fairly easy with a one-click deploy system in the 

First go the the AI platform page in the Google Cloud console.
https://console.cloud.google.com/ai-platform/notebooks/instances

Next you can deploy by hitting the 'New Instance' button and selecting 'Python'. There is no europe-west4 yet, so you can select europe-west2 (London) as you regions. For the machine type, leave it with the default setting so that we do not hit our computational quota.

Create an instance with you name and once its completed, hit the OPEN JUPYTERLAB button to connect to the instance.

We have to initialze the same google login command that we did for the prevous exercise, so open a new Terminal and see if you can login to your Google Cloud account, list the contents of the google cloud bucket where the manifest.json file is located and download the VCF file.







