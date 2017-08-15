# Overview

This project aims to address the difficulties in ingesting per-sample VCF files, provided for certain studies, into the [EVA-pipeline](https://github.com/EBIvariation/eva-pipeline/). 

As it currently stands, EVA-pipeline can only ingest "full" VCF files wherein information for all the samples is provided in a single VCF file. However, a single VCF file cannot be trivially generated for certain studies with a large number of samples, for ex: the [3000 Rice Genomes project](http://www.ebi.ac.uk/eva/?eva-study=PRJEB13618). In such cases, a solution, that can ingest these "single-sample" VCF files in parallel, is needed.

For a detailed discussion, see Proposal#2 in this [page](https://www.ebi.ac.uk/seqdb/confluence/display/VAR/Merging+Single-Sample+files).
      
# Project files

SingleSampleMerge.py - Apache Spark job that processes the individual study files in parallel, filters out monomorphic references and inserts the variants into Apache Cassandra.

ProcessVariantMatches.py - Apache Spark job that scans the individual study for the unique variant positions obtained from Cassandra (after SingleSampleMerge.py has run). Based on how many positions match in the individual sample files, a default genotype is assumed (ex:0/0 for >50% match, ./. otherwise) for the entire sample in order to minimize the number of inserts made into Cassandra.

spark_cluster_creation - Ansible playbook to create a Spark cluster.

# Steps to run the project

1. Create an Apache Spark Cluster as follows:
   1. Determine the master node. This is where you will download the project files and install Ansible. 
   2. Install Ansible on the master node with the instructions [here](http://docs.ansible.com/ansible/latest/intro_installation.html).
   3. Edit the **spark-hosts** file in the "spark_cluster_creation" folder with the IP addresses of the master node and the slave nodes.
   4. Run the Ansible playbook to deploy a Spark Cluster: ```ansible-playbook -i spark-hosts spark-install.yml```.

2. Create an Apache Cassandra cluster as described [here](https://github.com/EBIvariation/eva-variant-warehouse-research/blob/master/Cassandra_Evaluation/ansible_playbook/readme.md).
 
3. Download the study VCF files to a shared directory that can be accessed by all the Spark slave nodes.

4. Download the static binaries for bcftools from **/nfs/production3/eva/software** to a shared directory that can be accessed by all the Spark slave nodes.

5. Run SingleSampleMerge.py with arguments as follows:

   ```export SPARK_LOCAL_IP=<MASTER_NODE_IP> && ~/spark-2.2.0-bin-hadoop2.7/spark-submit --packages com.datastax.spark:spark-cassandra-connector_2.11:2.0.1 SingleSampleMerge.py <Study PRJ ID> <Study directory> <Cassandra Node IP1> <Cassandra Node IP2> <BCFTools directory>```
   
   Example:
   
   ```export SPARK_LOCAL_IP=192.168.0.14 && ~/spark-2.2.0-bin-hadoop2.7/bin/spark-submit --packages com.datastax.spark:spark-cassandra-connector_2.11:2.0.1 SingleSampleMerge.py PRJEB21300 /mnt/glusterVol/mergevcfinput 192.168.0.18 192.168.0.23 /mnt/glusterVol/bcftools```

6. If SingleSampleMerge.py ran successfully, run the ProcessVariantMatches.py as follows:

   ```export SPARK_LOCAL_IP=<MASTER_NODE_IP> && ~/spark-2.2.0-bin-hadoop2.7/spark-submit --packages com.datastax.spark:spark-cassandra-connector_2.11:2.0.1 ProcessVariantMatches.py <Study PRJ ID> <Default Genotype> <Alternate Genotype> <Full Path to study files> <Cassandra node IP1> <Cassandra node IP2> <BCF Tools Directory>```
   
   Example:
   
   ```export SPARK_LOCAL_IP=192.168.0.14 && ~/spark-2.2.0-bin-hadoop2.7/spark-submit --packages com.datastax.spark:spark-cassandra-connector_2.11:2.0.1 ProcessVariantMatches.py PRJEB21300 0/0 ./.  /mnt/glusterVol/mergevcfinput 192.168.0.18 192.168.0.23 /mnt/glusterVol/bcftools```
   
7. After the above steps are run, the following output tables are created in Cassandra in the "variant_ksp" keyspace:

    1. variant_ksp.variants_\<Study PRJ ID\> (ex: variant_ksp.variants_PRJEB21300)
    [variants]!(https://pasteboard.co/GFMiO2e.png)
    2. variant_ksp.headers_\<Study PRJ ID\> (ex: variant_ksp.headers_PRJEB21300)
    3. variant_ksp.sample_insert_log
    4. variant_ksp.study_info_\<Study PRJ ID\> (ex: variant_ksp.study_info_PRJEB21300)