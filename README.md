## GCP Workshop 

### Using the console

Your first entry point to GCP will likely be the console. It gives you a tool to both interact and monitor various GCP services and even
create a quick shell for more advances operations.

Hit the following link to start up the console: https://console.cloud.google.com/

### Projects and Accounts



### Using the command line

The fastest way to get to a shell is to use the "cloud shell" functionality. The cloud shell already has all tools installed and will 
default to the project you have open. It works by provisioning a small VM.

![Cloud Shell](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/cloud-shell.png)

Open a cloud shell and try the following:

```
gcloud auth login
gsutil ls
```

That said, the shell still runs in the browser so can be a little slow and not suited to long running interactions. You can also do 
everything from a terminal on your own maching with the `gcloud` and `gsutil` tools.

Open a terminal and try the following:
```
gcloud auth login
gsutil ls
```

### Google Cloud Storage (GCS)

Storage in the cloud is one of its most power tools and concepts. GCS is an object storage system which works with buckets, paths and files.

The first time, its best to create a bucket via the console:

![Create a bucket](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/storage-demo-1.png)

Important concepts:
* The region is where your bucket will be created. Best to pick a local region.
* The storage class impacts how your data is stored. Standard will be used for frequently accessed data. Nearline and Coldine for less 
accessed.
* The access control can be per bucket or per file, normally per bucket is fine.
* By default data in encrypted with Google's internal keys. Customer Managed Encryption is also available, but does carry a large 
performance penalty.

Create a bucket called `{yourname}-gcpdemo`. Bucket names have some restrictions: they must be globally unique; they must not have 
underscores; they must be lower-case. This is because you can expose data via direct URL access and must conform to the rules of DNS.

From there its best to move to the command line. After you've create your bucket, you can use the `gsutil` tool to upload, download, and
copy data between buckets. Try the following:

```
touch myfile.txt
gsutil cp myfile.txt gs://{yourname}-gcpdemo/
gsutil ls gs://{yourname}-gcpdemo/
```

You can also browse your data from within the console. Find Storage in your console and you can see your new bucket and file there.

You can share your data with other users using the Access Control List of the bucket. Adding users is easiest via the console, but can also 
be done with `gsutil`. Give it a try with your neighbour!

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

You can also create and manage vms via the command line. 

```bash
gcloud instances create my-instance --zone=europe-west4-a --boot-disk-size 100
```

The [VM instance overview](https://console.cloud.google.com/compute/instances) in the console is a nice way to see the status of all VMs:

![VM Instances](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/compute-demo-2.png)

The SSH button gives you a quick browser shell to the machine. Again, the browser has its limitations, so sometimes it's better to tunnel 
into your VM from your own terminal:

```bash
gcloud beta compute --project "your-project" ssh --zone "europe-west4-a" "your-vm"
```

Try SSH'ing into your new VM and running the same `gsutil` commands we ran in the previous sections.

### 
 







