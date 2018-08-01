# RNA-Seq Workflow Tutorial

### 1. Make project directory
```
# the project directory contains specific GROUP name and current DATE
GROUP=XXX
DATE=`date +"%Y%m%d"`
mkdir ~/Project/${GROUP}_${DATE}
```

### 2. Clone the repository
```
cd ~/Project/${GROUP}_${DATE}
git clone https://github.com/bioxfu/RNA-Seq
cd RNA-Seq
```

### 3. Copy/Download the raw data
```
DATAPATH=/the/path/of/the/raw/data/on/HPC

# if you are working on the HPC, copy the raw data
./script/copy_rawdata.sh $DATAPATH

# if you are working on the local machine, download the raw data
./script/copy_rawdata.sh $DATAPATH --download
```

### 4. Rename the raw data
```
# dry run to check if mv command is correct
./script/rename_rawdata.sh --dry-run

# then do it 
./script/rename_rawdata.sh
```

### 5. Create *config.yaml* and *Snakefile* based on the examples
```
cp example/example.HPC.plant.config.yaml config.yaml
cp example/example.plant.Snakefile Snakefile
cp example/example_design_table.tsv design_table.tsv

# edit config.yaml and design_table.tsv
```

### 6. Initiate the project
```
source init.sh
```

### 7. Dry run the workflow to check any mistakes
```
./dry_run.sh
```

### 8. Start the workflow
```
# if you are working on the HPC
./run_HPC.sh

# if you are working on the local machine
./run.sh

# check the workflow progress in nohup.out file
tail nohup.log 
```

### 9. Remove the temporary files
```
./clean.sh
```

