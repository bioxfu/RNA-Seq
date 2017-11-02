# RNA-Seq workflow

### 1. Clone the repository
```
git clone https://github.com/bioxfu/RNA-Seq
```

### 2. Create a new branch (*project*)
```
git branch project
```

### 3. Checkout *project* branch
```
git checkout project
```

### 4. Create *config.yaml* and *Snakefile* based on the examples

### 5. Initiate the project
```
source init.sh
```

### 6. Dry run the workflow to check any mistakes
```
./dry_run.sh
```

### 7. Start the workflow
```
./run.sh
```

### 8. Check the workflow progress in *nohup.out* 

### 9. Remove the temporary files
```
./clean.sh
```

