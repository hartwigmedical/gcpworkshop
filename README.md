## GCP Workshop 

### Using the console

Your first entry point to GCP will likely be the console. It gives you a tool to both interact and monitor various GCP services and even
create a quick shell for more advances operations.

Hit the following link to start up the console: https://console.cloud.google.com/

### Using the command line

The fastest way to get to a shell is to use the "cloud shell" functionality. The cloud shell already has all tools installed and will 
default to the project you have open. It works by provisioning a small VM.

TODO: CS screenshot

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

![Create a bucket](https://github.com/hartwigmedical/gcpworkshop/images/)

Important concepts:
* The region is where your bucket will be created. Best to pick a local region.
* The storage class impacts how your data is stored. 






