# AIAP Test Version

## Usage

### 0. Install Nextflow
Refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/install.html) for installation instructions.

### 1. Build the Singularity Image
First, build the singularity image using the provided definition file:
```
singularity build aiap-nf.sif aiap-nf.def
```
Alternatively, you can use any other singularity image that contains the required tools.

### 2. Prepare Input Files
- **Sample Sheet**: Edit `samplesheet-demo.csv` (filename can be changed) to include your samples.
- **Configuration**: Edit `params-demo.json` (filename can be changed) to specify your pipeline parameters.

### 3. Run the Pipeline
You can run the pipeline using:
```
nextflow run tyfu-cc/aiap -params-file /path/to/params/json/file -with-singularity /path/to/your/sif/file
```
Or, if you've cloned the pipeline locally, navigate to its directory and run:
```
nextflow run main.nf -params-file /path/to/your/params/json/file -with-singularity /path/to/your/sif/file
```
