## run minigraph cactus for selected citrus genomes
https://github.com/ComparativeGenomicsToolkit/cactus/tree/master
https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md
minigraph cactus for aligning the same species
Progressive Cactus for aligning different species

## install cactus 
https://github.com/ComparativeGenomicsToolkit/cactus/releases
```
singularity build cactus29.sif docker:quay.io/comparative-genomics-toolkit/cactus:v2.9.0
```

## test run cactus
https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/mc-pangenomes/README.md
### location: 
(base) yongjia@ChengdaoServer:/media/yongjia/Elements/yongjia/cactus_test$
### download 4-t2t-orangs-mc-2023v2.seqfile
```
mPonAbe1_pri	https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/655/GCA_028885655.2/GCA_028885655.2.fa.gz
mPonAbe1_alt	https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/685/GCA_028885685.2/GCA_028885685.2.fa.gz
mPonPyg2.1	https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/625/GCA_028885625.2/GCA_028885625.2.fa.gz
mPonPyg2.2	https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/525/GCA_028885525.2/GCA_028885525.2.fa.gz
```
### command cactus.sh
```
cactus-pangenome ./js-pg ./4-t2t-orangs-mc-2023v2.seqfile --outDir 4-t2t-orangs-mc-2023v2 --outName 4-t2t-orangs-mc-2023v2 --reference mPonAbe1_pri mPonAbe1_alt \
--noSplit --gbz clip full --gfa clip full --xg clip full --odgi --vcf --giraffe clip --haplo clip --vcfReference mPonAbe1_pri mPonAbe1_alt --logFile 4-t2t-orangs-mc-2023v2.log  \
--coordinationDir /data/tmp --batchLogsDir ./batch-logs --consMemory 230Gi --indexMemory 230Gi --mgMemory 230Gi --mgCores 60 --mapCores 8 --consCores 60 --indexCores 60 --giraffe clip
```
```
cactus-pangenome
./js-pg
./4-t2t-orangs-mc-2023v2.seqfile
--outDir
4-t2t-orangs-mc-2023v2
--outName
4-t2t-orangs-mc-2023v2
--reference
mPonAbe1_pri
mPonAbe1_alt
--noSplit
--gbz
clip
full
--gfa
clip
full
--xg
clip
full
--odgi
--vcf
--giraffe
clip
--haplo
clip
--vcfReference
mPonAbe1_pri
mPonAbe1_alt
--logFile
4-t2t-orangs-mc-2023v2.log
--coordinationDir
/data/tmp
--batchLogsDir
./batch-logs
--consMemory
230Gi
--indexMemory
230Gi
--mgMemory
230Gi
--mgCores
60
--mapCores
8
--consCores
60
--indexCores
60
--giraffe
clip
```
### run using singularity in tmux
```
singularity shell -B ${PWD}:/data /media/yongjia/Elements/yongjia/containers/cactus29.sif
bash cactus.sh &> log.txt
```
## Notes
```
```

## cactus-pangenome usage
```
Singularity> cactus-pangenome -h
usage: cactus-pangenome [-h] [--config PATH] [--logCritical] [--logError] [--logWarning] [--logDebug] [--logInfo] [--logOff] [--logLevel {Critical,Error,Warning,Debug,Info,critical,error,warning,debug,info,CRITICAL,ERROR,WARNING,DEBUG,INFO}] [--logFile LOGFILE]
                        [--rotatingLogging] [--workDir PATH] [--coordinationDir PATH] [--noStdOutErr] [--stats] [--clean {always,onError,never,onSuccess}] [--cleanWorkDir {always,onError,never,onSuccess}] [--clusterStats [OPT_PATH]] [--restart]
                        [--batchSystem {aws_batch,single_machine,grid_engine,lsf,mesos,slurm,torque,htcondor,kubernetes}] [--disableHotDeployment] [--disableAutoDeployment DISABLEAUTODEPLOYMENT] [--maxJobs MAX_JOBS] [--maxLocalJobs MAX_LOCAL_JOBS] [--manualMemArgs]
                        [--runLocalJobsOnWorkers] [--coalesceStatusCalls] [--statePollingWait STATEPOLLINGWAIT] [--batchLogsDir BATCH_LOGS_DIR] [--awsBatchRegion AWS_BATCH_REGION] [--awsBatchQueue AWS_BATCH_QUEUE] [--awsBatchJobRoleArn AWS_BATCH_JOB_ROLE_ARN]
                        [--scale SCALE] [--dont_allocate_mem | --allocate_mem] [--symlinkImports BOOL] [--moveOutputs BOOL] [--caching BOOL] [--provisioner {aws,gce,None}] [--nodeTypes NODETYPES] [--maxNodes INT[,INT...]] [--minNodes INT[,INT...]] [--targetTime INT]
                        [--betaInertia FLOAT] [--scaleInterval INT] [--preemptibleCompensation FLOAT] [--nodeStorage INT] [--nodeStorageOverrides NODETYPE:NODESTORAGE[,NODETYPE:NODESTORAGE...]] [--metrics BOOL] [--assumeZeroOverhead BOOL] [--maxServiceJobs INT]
                        [--maxPreemptibleServiceJobs INT] [--deadlockWait INT] [--deadlockCheckInterval INT] [--defaultMemory DEFAULTMEMORY] [--defaultCores FLOAT] [--defaultDisk INT] [--defaultAccelerators ACCELERATOR[,ACCELERATOR...]] [--defaultPreemptible [BOOL]]
                        [--maxCores INT] [--maxMemory INT] [--maxDisk INT] [--retryCount INT] [--enableUnlimitedPreemptibleRetries BOOL] [--doubleMem BOOL] [--maxJobDuration INT] [--rescueJobsFrequency INT] [--maxLogFileSize MAXLOGFILESIZE] [--writeLogs [OPT_PATH]]
                        [--writeLogsGzip [OPT_PATH]] [--writeLogsFromAllJobs BOOL] [--writeMessages PATH] [--realTimeLogging REALTIMELOGGING] [--disableChaining BOOL] [--disableJobStoreChecksumVerification BOOL] [--sseKey PATH] [--setEnv NAME=VALUE or NAME]
                        [--servicePollingInterval FLOAT] [--forceDockerAppliance BOOL] [--statusWait INT] [--disableProgress BOOL] [--debugWorker] [--disableWorkerOutputCapture] [--badWorker FLOAT] [--badWorkerFailInterval FLOAT] --outDir OUTDIR --outName OUTNAME
                        --reference REFERENCE [REFERENCE ...] [--mgCores MGCORES] [--mgMemory MGMEMORY] [--mapCores MAPCORES] [--delFilter DELFILTER] [--refContigs [REFCONTIGS ...]] [--otherContig OTHERCONTIG] [--permissiveContigFilter [PERMISSIVECONTIGFILTER]]
                        [--noSplit] [--maxLen MAXLEN] [--consCores CONSCORES] [--consMemory CONSMEMORY] [--clip CLIP] [--filter FILTER] [--gfa [GFA ...]] [--gbz [GBZ ...]] [--xg [XG ...]] [--odgi [ODGI ...]] [--viz [VIZ ...]] [--draw [DRAW ...]]
                        [--chrom-vg [CHROM_VG ...]] [--chrom-og [CHROM_OG ...]] [--vcf [VCF ...]] [--vcfReference VCFREFERENCE [VCFREFERENCE ...]] [--vcfbub VCFBUB] [--giraffe [GIRAFFE ...]] [--haplo [HAPLO ...]] [--indexCores INDEXCORES] [--indexMemory INDEXMEMORY]
                        [--chop [CHOP]] [--configFile CONFIGFILE] [--latest] [--containerImage CONTAINERIMAGE] [--binariesMode {docker,local,singularity}]
                        jobStore seqFile

positional arguments:
  seqFile               Seq file (will be modified if necessary to include graph Fasta sequence)

options:
  -h, --help            show this help message and exit
  --outDir OUTDIR       Output directory
  --outName OUTNAME     Output name (without extension)
  --reference REFERENCE [REFERENCE ...]
                        Reference event name(s). The first will be the "true" reference and will be left unclipped and uncollapsed. It also should have been used with --reference in all upstream commands. Other names will be promoted to reference paths in vg
  --mgCores MGCORES     Number of cores for minigraph construction (defaults to the same as --maxCores).
  --mgMemory MGMEMORY   Memory in bytes for the minigraph construction job (defaults to an estimate based on the input data size). Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))
  --mapCores MAPCORES   Number of cores for minigraph map. Overrides graphmap cpu in configuration
  --delFilter DELFILTER
                        Filter out split-mapping-implied deletions > Nbp (default will be "delFilter" from the config
  --refContigs [REFCONTIGS ...]
                        Subset to these reference contigs (multiple allowed)
  --otherContig OTHERCONTIG
                        Lump all reference contigs unselected by above options into single one with this name
  --permissiveContigFilter [PERMISSIVECONTIGFILTER]
                        If specified, override the configuration to accept contigs so long as they have at least given fraction of coverage (0.25 if no fraction specified). This can increase sensitivity of very small, fragmented and/or diverse assemblies.
  --noSplit             Do not split by ref chromsome. This will require much more memory and potentially produce a more complex graph
  --maxLen MAXLEN       Only align up to this many bases (overrides <bar bandingLimit> and <caf maxRecoverableChainLength> in configuration)[default=10000]
  --consCores CONSCORES
                        Number of cores for each cactus_consolidated job (defaults to all available / maxCores on single_machine)
  --consMemory CONSMEMORY
                        Memory in bytes for each cactus_consolidated job (defaults to an estimate based on the input data size). Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))
  --clip CLIP           Generate clipped graph by removing anything longer than this amount that is unaligned to the underlying minigraph. Set to 0 to disable (must also set --filter 0 as well). [default=10000]
  --filter FILTER       Generate a frequency filtered graph (from the clipped graph) by removing any sequence present in fewer than this many sequences. Set to 0 to disable. [default=2]
  --gfa [GFA ...]       Produce a GFA for given graph type(s) if specified. Valid types are 'full', 'clip', and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space. [--gfa clip
                        assumed by default]
  --gbz [GBZ ...]       Generate GBZ/snarls indexes for the given graph type(s) if specified. Valid types are 'full', 'clip' and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space.
                        --giraffe will also produce these (and other) indexes
  --xg [XG ...]         Generate XG index (from GBZ) for the given graph type(s) if specified. Valid types are 'full', 'clip' and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space.
  --odgi [ODGI ...]     Generate ODGI (.og) versions for the given graph type(s) if specified. Valid types are 'full', 'clip' and 'filter'. If no type specified 'full' will be used. Multiple types can be provided separated by a space.
  --viz [VIZ ...]       Generate 1D ODGI visualizations of each chromosomal graph for the given graph type(s) if specified. Valid types are 'full', 'clip' (but NOT `filter`). If no type specified 'full' will be used. Multiple types can be provided separated by a space.
  --draw [DRAW ...]     Generate 2D ODGI visualizations of each chromosomal graph for the given graph type(s) if specified. WARNING: EXPERIMENTAL OPTION. More work needs to be done to figure out the best ODGI parameters to use, and it can be quite slow in the meantime.
                        For larger graphs, use --chrom-og and run the drawing by hand in order to avoid having draw issues prevent you from getting the rest of your output. Valid types are 'full', 'clip' (but NOT `filter`). If no type specified 'full' will be used.
                        Multiple types can be provided separated by a space.
  --chrom-vg [CHROM_VG ...]
                        Produce a directory of chromosomal graphs is vg format for the graph type(s) specified. Valid typs are 'full', 'clip' and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided
                        separated by a space. The output will be the <outDir>/<outName>.chroms/ directory
  --chrom-og [CHROM_OG ...]
                        Produce a directory of chromosomal graphs is odgi format for the graph type(s) specified. Valid typs are 'full', 'clip' and 'filter'. If no type specified 'full' will be used. Multiple types can be provided separated by a space. The output will
                        be the <outDir>/<outName>.chroms/ directory
  --vcf [VCF ...]       Generate a VCF from the given graph type(s). Valid types are 'full', 'clip' and 'filter'. If no type specified, 'clip' will be used ('full' used if clipping disabled). Multipe types can be provided separated by space
  --vcfReference VCFREFERENCE [VCFREFERENCE ...]
                        If multiple references were provided with --reference, this option can be used to specify a subset for vcf creation with --vcf. By default, --vcf will create VCFs for the first reference only
  --vcfbub VCFBUB       Use vcfbub to flatten nested sites (sites with reference alleles > this will be replaced by their children)). Setting to 0 will disable, only prudcing full VCF [default=100000].
  --giraffe [GIRAFFE ...]
                        Generate Giraffe (.dist, .min) indexes for the given graph type(s). Valid types are 'full', 'clip' and 'filter'. If not type specified, 'filter' will be used (will fall back to 'clip' than full if filtering, clipping disabled, respectively).
                        Multiple types can be provided seperated by a space
  --haplo [HAPLO ...]   Generate haplotype subsampling (.ri, .hapl) indexes for the given graph type(s). Haplotype subsampling is a new, better alternative filtering by allele frequency. Valid types are 'full' and 'clip'. If not type specified, 'clip' will be used
                        ('full' will be used if clipping disabled). Multiple types can be provided seperated by a space. Must be used in conjunction with --giraffe
  --indexCores INDEXCORES
                        cores for general indexing and VCF constructions (defaults to the same as --maxCores)
  --indexMemory INDEXMEMORY
                        Memory in bytes for each indexing and vcf construction job job (defaults to an estimate based on the input data size). If specified will also be used to upper-bound per-chromosome memory estimates -- ie no job will request more than this much
                        memory.Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))
  --chop [CHOP]         chop all graph nodes to be at most this long (default=1024 specified without value). By default, nodes are only chopped for GBZ-derived formats, but left unchopped in the GFA, VCF, etc. If this option is used, the GBZ and GFA should have
                        consistent IDs
  --configFile CONFIGFILE
                        Specify cactus configuration file
  --latest              Use the latest version of the docker container rather than pulling one matching this version of cactus
  --containerImage CONTAINERIMAGE
                        Use the the specified pre-built containter image rather than pulling one from quay.io
  --binariesMode {docker,local,singularity}
                        The way to run the Cactus binaries

  --config PATH         Get options from a config file.

Logging Options:
  --logCritical         Turn on loglevel Critical. Default: INFO.
  --logError            Turn on loglevel Error. Default: INFO.
  --logWarning          Turn on loglevel Warning. Default: INFO.
  --logDebug            Turn on loglevel Debug. Default: INFO.
  --logInfo             Turn on loglevel Info. Default: INFO.
  --logOff              Same as --logCRITICAL.
  --logLevel {Critical,Error,Warning,Debug,Info,critical,error,warning,debug,info,CRITICAL,ERROR,WARNING,DEBUG,INFO}
                        Set the log level. Default: INFO. Options: ['Critical', 'Error', 'Warning', 'Debug', 'Info', 'critical', 'error', 'warning', 'debug', 'info', 'CRITICAL', 'ERROR', 'WARNING', 'DEBUG', 'INFO'].
  --logFile LOGFILE     File to log in.
  --rotatingLogging     Turn on rotating logging, which prevents log files from getting too big.

Toil core options.:
  Options to specify the location of the Toil workflow and turn on stats collation about the performance of jobs.

  jobStore              The location of the job store for the workflow. A job store holds persistent information about the jobs, stats, and files in a workflow. If the workflow is run with a distributed batch system, the job store must be accessible by all worker nodes.
                        Depending on the desired job store implementation, the location should be formatted according to one of the following schemes: file:<path> where <path> points to a directory on the file systen aws:<region>:<prefix> where <region> is the name of
                        an AWS region like us-west-2 and <prefix> will be prepended to the names of any top-level AWS resources in use by job store, e.g. S3 buckets. google:<project_id>:<prefix> TODO: explain For backwards compatibility, you may also specify ./foo
                        (equivalent to file:./foo or just file:foo) or /bar (equivalent to file:/bar).
  --workDir PATH        Absolute path to directory where temporary files generated during the Toil run should be placed. Standard output and error from batch system jobs (unless --noStdOutErr is set) will be placed in this directory. A cache directory may be placed in
                        this directory. Temp files and folders will be placed in a directory toil-<workflowID> within workDir. The workflowID is generated by Toil and will be reported in the workflow logs. Default is determined by the variables (TMPDIR, TEMP, TMP) via
                        mkdtemp. This directory needs to exist on all machines running jobs; if capturing standard output and error from batch system jobs is desired, it will generally need to be on a shared file system. When sharing a cache between containers on a
                        host, this directory must be shared between the containers.
  --coordinationDir PATH
                        Absolute path to directory where Toil will keep state and lock files.When sharing a cache between containers on a host, this directory must be shared between the containers.
  --noStdOutErr         Do not capture standard output and error from batch system jobs.
  --stats               Records statistics about the toil workflow to be used by 'toil stats'.
  --clean {always,onError,never,onSuccess}
                        Determines the deletion of the jobStore upon completion of the program. Choices: ['always', 'onError', 'never', 'onSuccess']. The --stats option requires information from the jobStore upon completion so the jobStore will never be deleted with
                        that flag. If you wish to be able to restart the run, choose 'never' or 'onSuccess'. Default is 'never' if stats is enabled, and 'onSuccess' otherwise.
  --cleanWorkDir {always,onError,never,onSuccess}
                        Determines deletion of temporary worker directory upon completion of a job. Choices: ['always', 'onError', 'never', 'onSuccess']. Default = always. WARNING: This option should be changed for debugging only. Running a full pipeline with this
                        option could fill your disk with excessive intermediate data.
  --clusterStats [OPT_PATH]
                        If enabled, writes out JSON resource usage statistics to a file. The default location for this file is the current working directory, but an absolute path can also be passed to specify where this file should be written. This options only applies
                        when using scalable batch systems.

Toil options for restarting an existing workflow.:
  Allows the restart of an existing workflow

  --restart             If --restart is specified then will attempt to restart existing workflow at the location pointed to by the --jobStore option. Will raise an exception if the workflow does not exist

Toil options for specifying the batch system.:
  Allows the specification of the batch system.

  --batchSystem {aws_batch,single_machine,grid_engine,lsf,mesos,slurm,torque,htcondor,kubernetes}
                        The type of batch system to run the job(s) with, currently can be one of aws_batch, single_machine, grid_engine, lsf, mesos, slurm, torque, htcondor, kubernetes. Default = single_machine
  --disableHotDeployment
                        Hot-deployment was renamed to auto-deployment. Option now redirects to --disableAutoDeployment. Left in for backwards compatibility.
  --disableAutoDeployment DISABLEAUTODEPLOYMENT
                        Should auto-deployment of Toil Python workflows be deactivated? If True, the workflow's Python code should be present at the same location on all workers. Default = False.
  --maxJobs MAX_JOBS    Specifies the maximum number of jobs to submit to the backing scheduler at once. Not supported on Mesos or AWS Batch. Use 0 for unlimited. Defaults to unlimited.
  --maxLocalJobs MAX_LOCAL_JOBS
                        Specifies the maximum number of housekeeping jobs to run sumultaneously on the local system. Use 0 for unlimited. Defaults to the number of local cores (64).
  --manualMemArgs       Do not add the default arguments: 'hv=MEMORY' & 'h_vmem=MEMORY' to the qsub call, and instead rely on TOIL_GRIDGENGINE_ARGS to supply alternative arguments. Requires that TOIL_GRIDGENGINE_ARGS be set.
  --runLocalJobsOnWorkers, --runCwlInternalJobsOnWorkers
                        Whether to run jobs marked as local (e.g. CWLScatter) on the worker nodes instead of the leader node. If false (default), then all such jobs are run on the leader node. Setting this to true can speed up CWL pipelines for very large workflows with
                        many sub-workflows and/or scatters, provided that the worker pool is large enough.
  --coalesceStatusCalls
                        Ask for job statuses from the batch system in a batch. Deprecated; now always enabled where supported.
  --statePollingWait STATEPOLLINGWAIT
                        Time, in seconds, to wait before doing a scheduler query for job state. Return cached results if within the waiting period. Only works for grid engine batch systems such as gridengine, htcondor, torque, slurm, and lsf.
  --batchLogsDir BATCH_LOGS_DIR
                        Directory to tell the backing batch system to log into. Should be available on both the leader and the workers, if the backing batch system writes logs to the worker machines' filesystems, as many HPC schedulers do. If unset, the Toil work
                        directory will be used. Only works for grid engine batch systems such as gridengine, htcondor, torque, slurm, and lsf.
  --awsBatchRegion AWS_BATCH_REGION
                        The AWS region containing the AWS Batch queue to submit to.
  --awsBatchQueue AWS_BATCH_QUEUE
                        The name or ARN of the AWS Batch queue to submit to.
  --awsBatchJobRoleArn AWS_BATCH_JOB_ROLE_ARN
                        The ARN of an IAM role to run AWS Batch jobs as, so they can e.g. access a job store. Must be assumable by ecs-tasks.amazonaws.com.
  --scale SCALE         A scaling factor to change the value of all submitted tasks's submitted cores. Used in the single_machine batch system. Useful for running workflows on smaller machines than they were designed for, by setting a value less than 1. (default: 1)
  --dont_allocate_mem   A flag that can block allocating memory with '--mem' for job submissions on SLURM since some system servers may reject any job request that explicitly specifies the memory allocation. The default is to always allocate memory.
  --allocate_mem        A flag that can block allocating memory with '--mem' for job submissions on SLURM since some system servers may reject any job request that explicitly specifies the memory allocation. The default is to always allocate memory.

Toil options for configuring storage.:
  Allows configuring Toil's data storage.

  --symlinkImports BOOL
                        When using a filesystem based job store, CWL input files are by default symlinked in. Setting this option to True instead copies the files into the job store, which may protect them from being modified externally. When set to False, as long as
                        caching is enabled, Toil will protect the file automatically by changing the permissions to read-only.default=True
  --moveOutputs BOOL    When using a filesystem based job store, output files are by default moved to the output directory, and a symlink to the moved exported file is created at the initial location. Setting this option to True instead copies the files into the output
                        directory. Applies to filesystem-based job stores only.default=False
  --caching BOOL        Enable or disable caching for your workflow, specifying this overrides default from job store

Toil options for autoscaling the cluster of worker nodes.:
  Allows the specification of the minimum and maximum number of nodes in an autoscaled cluster, as well as parameters to control the level of provisioning.

  --provisioner {aws,gce,None}, -p {aws,gce,None}
                        The provisioner for cluster auto-scaling. This is the main Toil '--provisioner' option, and defaults to None for running on single machine and non-auto-scaling batch systems. The currently supported choices are ['aws', 'gce', None]. The default
                        is None.
  --nodeTypes NODETYPES
                        Specifies a list of comma-separated node types, each of which is composed of slash-separated instance types, and an optional spot bid set off by a colon, making the node type preemptible. Instance types may appear in multiple node types, and the
                        same node type may appear as both preemptible and non-preemptible. Valid argument specifying two node types: c5.4xlarge/c5a.4xlarge:0.42,t2.large Node types: c5.4xlarge/c5a.4xlarge:0.42 and t2.large Instance types: c5.4xlarge, c5a.4xlarge, and
                        t2.large Semantics: Bid $0.42/hour for either c5.4xlarge or c5a.4xlarge instances, treated interchangeably, while they are available at that price, and buy t2.large instances at full price. default=[]
  --maxNodes INT[,INT...]
                        Maximum number of nodes of each type in the cluster, if using autoscaling, provided as a comma-separated list. The first value is used as a default if the list length is less than the number of nodeTypes. default=[10]
  --minNodes INT[,INT...]
                        Mininum number of nodes of each type in the cluster, if using auto-scaling. This should be provided as a comma-separated list of the same length as the list of node types. default=[0]
  --targetTime INT      Sets how rapidly you aim to complete jobs in seconds. Shorter times mean more aggressive parallelization. The autoscaler attempts to scale up/down so that it expects all queued jobs will complete within targetTime seconds. default=1800
  --betaInertia FLOAT   A smoothing parameter to prevent unnecessary oscillations in the number of provisioned nodes. This controls an exponentially weighted moving average of the estimated number of nodes. A value of 0.0 disables any smoothing, and a value of 0.9 will
                        smooth so much that few changes will ever be made. Must be between 0.0 and 0.9. default=0.1
  --scaleInterval INT   The interval (seconds) between assessing if the scale of the cluster needs to change. default=60
  --preemptibleCompensation FLOAT, --preemptableCompensation FLOAT
                        The preference of the autoscaler to replace preemptible nodes with non-preemptible nodes, when preemptible nodes cannot be started for some reason. This value must be between 0.0 and 1.0, inclusive. A value of 0.0 disables such compensation, a
                        value of 0.5 compensates two missing preemptible nodes with a non-preemptible one. A value of 1.0 replaces every missing pre-emptable node with a non-preemptible one. default=0.0
  --nodeStorage INT     Specify the size of the root volume of worker nodes when they are launched in gigabytes. You may want to set this if your jobs require a lot of disk space. (default=50).
  --nodeStorageOverrides NODETYPE:NODESTORAGE[,NODETYPE:NODESTORAGE...]
                        Comma-separated list of nodeType:nodeStorage that are used to override the default value from --nodeStorage for the specified nodeType(s). This is useful for heterogeneous jobs where some tasks require much more disk than others.
  --metrics BOOL        Enable the prometheus/grafana dashboard for monitoring CPU/RAM usage, queue size, and issued jobs.
  --assumeZeroOverhead BOOL
                        Ignore scheduler and OS overhead and assume jobs can use every last byte of memory and disk on a node when autoscaling.

Toil options for limiting the number of service jobs and detecting service deadlocks:
  Allows the specification of the maximum number of service jobs in a cluster. By keeping this limited we can avoid nodes occupied with services causing deadlocks.

  --maxServiceJobs INT  The maximum number of service jobs that can be run concurrently, excluding service jobs running on preemptible nodes. default=9223372036854775807
  --maxPreemptibleServiceJobs INT
                        The maximum number of service jobs that can run concurrently on preemptible nodes. default=9223372036854775807
  --deadlockWait INT    Time, in seconds, to tolerate the workflow running only the same service jobs, with no jobs to use them, before declaring the workflow to be deadlocked and stopping. default=60
  --deadlockCheckInterval INT
                        Time, in seconds, to wait between checks to see if the workflow is stuck running only service jobs, with no jobs to use them. Should be shorter than --deadlockWait. May need to be increased if the batch system cannot enumerate running jobs
                        quickly enough, or if polling for running jobs is placing an unacceptable load on a shared cluster.default=30

Toil options for cores/memory requirements.:
  The options to specify default cores/memory requirements (if not specified by the jobs themselves), and to limit the total amount of memory/cores requested from the batch system.

  --defaultMemory DEFAULTMEMORY
                        The default amount of memory to request for a job. Only applicable to jobs that do not specify an explicit value for this requirement. Standard suffixes like K, Ki, M, Mi, G or Gi are supported. Default is 2.0 Gi.
  --defaultCores FLOAT  The default amount of cpu to request for a job. Only applicable to jobs that do not specify an explicit value for this requirement. Fractions of a core (for example 0.1) are supported on some batch systems [mesos, single_machine]. Default is 1.
  --defaultDisk INT     The default amount of disk to request for a job. Only applicable to jobs that do not specify an explicit value for this requirement. Standard suffixes like K, Ki, M, Mi, G or Gi are supported. Default is 2.0 Gi.
  --defaultAccelerators ACCELERATOR[,ACCELERATOR...]
                        The default amount of accelerators to request for a job. Only applicable to jobs that do not specify an explicit value for this requirement. Each accelerator specification can have a type (gpu [default], nvidia, amd, cuda, rocm, opencl, or a
                        specific model like nvidia-tesla-k80), and a count [default: 1]. If both a type and a count are used, they must be separated by a colon. If multiple types of accelerators are used, the specifications are separated by commas. Default is [].
  --defaultPreemptible [BOOL], --defaultPreemptable [BOOL]
                        Make all jobs able to run on preemptible (spot) nodes by default.
  --maxCores INT        The max amount of cpu to request for a job. Only applicable to jobs that do not specify an explicit value for this requirement. Fractions of a core (for example 0.1) are supported on some batch systems [mesos, single_machine]. Default is
                        9223372036854775807.
  --maxMemory INT       The max amount of memory to request for a job. Only applicable to jobs that do not specify an explicit value for this requirement. Standard suffixes like K, Ki, M, Mi, G or Gi are supported. Default is 8.0 Ei.
  --maxDisk INT         The max amount of disk to request for a job. Only applicable to jobs that do not specify an explicit value for this requirement. Standard suffixes like K, Ki, M, Mi, G or Gi are supported. Default is 8.0 Ei.

Toil options for rescuing/killing/restarting jobs.:
  The options for jobs that either run too long/fail or get lost (some batch systems have issues!).

  --retryCount INT      Number of times to retry a failing job before giving up and labeling job failed. default=1
  --enableUnlimitedPreemptibleRetries BOOL, --enableUnlimitedPreemptableRetries BOOL
                        If set, preemptible failures (or any failure due to an instance getting unexpectedly terminated) will not count towards job failures and --retryCount.
  --doubleMem BOOL      If set, batch jobs which die to reaching memory limit on batch schedulers will have their memory doubled and they will be retried. The remaining retry count will be reduced by 1. Currently supported by LSF.
  --maxJobDuration INT  Maximum runtime of a job (in seconds) before we kill it (this is a lower bound, and the actual time before killing the job may be longer). default=9223372036854775807
  --rescueJobsFrequency INT
                        Period of time to wait (in seconds) between checking for missing/overlong jobs, that is jobs which get lost by the batch system. Expert parameter. default=60

Toil log management options.:
  Options for how Toil should manage its logs.

  --maxLogFileSize MAXLOGFILESIZE
                        The maximum size of a job log file to keep (in bytes), log files larger than this will be truncated to the last X bytes. Setting this option to zero will prevent any truncation. Setting this option to a negative value will truncate from the
                        beginning. Default=62.5 Ki
  --writeLogs [OPT_PATH]
                        Write worker logs received by the leader into their own files at the specified path. Any non-empty standard output and error from failed batch system jobs will also be written into files at this path. The current working directory will be used if
                        a path is not specified explicitly. Note: By default only the logs of failed jobs are returned to leader. Set log level to 'debug' or enable '--writeLogsFromAllJobs' to get logs back from successful jobs, and adjust 'maxLogFileSize' to control
                        the truncation limit for worker logs.
  --writeLogsGzip [OPT_PATH]
                        Identical to --writeLogs except the logs files are gzipped on the leader.
  --writeLogsFromAllJobs BOOL
                        Whether to write logs from all jobs (including the successful ones) without necessarily setting the log level to 'debug'. Ensure that either --writeLogs or --writeLogsGzip is set if enabling this option.
  --writeMessages PATH  File to send messages from the leader's message bus to.
  --realTimeLogging REALTIMELOGGING
                        Enable real-time logging from workers to leader

Toil miscellaneous options.:
  Everything else.

  --disableChaining BOOL
                        Disables chaining of jobs (chaining uses one job's resource allocation for its successor job if possible).
  --disableJobStoreChecksumVerification BOOL
                        Disables checksum verification for files transferred to/from the job store. Checksum verification is a safety check to ensure the data is not corrupted during transfer. Currently only supported for non-streaming AWS files.
  --sseKey PATH         Path to file containing 32 character key to be used for server-side encryption on awsJobStore or googleJobStore. SSE will not be used if this flag is not passed.
  --setEnv NAME=VALUE or NAME, -e NAME=VALUE or NAME
                        Set an environment variable early on in the worker. If VALUE is null, it will be looked up in the current environment. Independently of this option, the worker will try to emulate the leader's environment before running a job, except for some
                        variables known to vary across systems. Using this option, a variable can be injected into the worker process itself before it is started.
  --servicePollingInterval FLOAT
                        Interval of time service jobs wait between polling for the existence of the keep-alive flag. Default: 60.0
  --forceDockerAppliance BOOL
                        Disables sanity checking the existence of the docker image specified by TOIL_APPLIANCE_SELF, which Toil uses to provision mesos for autoscaling.
  --statusWait INT      Seconds to wait between reports of running jobs.
  --disableProgress BOOL
                        Disables the progress bar shown when standard error is a terminal.

Toil debug options.:
  Debug options for finding problems or helping with testing.

  --debugWorker         Experimental no forking mode for local debugging. Specifically, workers are not forked and stderr/stdout are not redirected to the log.
  --disableWorkerOutputCapture
                        Let worker output go to worker's standard out/error instead of per-job logs.
  --badWorker FLOAT     For testing purposes randomly kill --badWorker proportion of jobs using SIGKILL. default=0.0
  --badWorkerFailInterval FLOAT
                        When killing the job pick uniformly within the interval from 0.0 to --badWorkerFailInterval seconds after the worker starts. default=0.01
```
