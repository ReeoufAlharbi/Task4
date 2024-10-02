# Task4

### Problem 1

```shell
python Genome1.py Genome_File.fasta > output1_genome.txt
```

### Problem 2

```shell
python Genome2.py Genome_File.fasta > output2_genome.txt
```

### Problem 3

```shell
python Genome3.py Genome3.fasta > output3_genome.txt
```

### Problem 4 

```shell
Get-ChildItem *.fna | ForEach-Object {
>>     $outputFile = "$($_.BaseName)_output4.fna"
>>     $output = python Genome4.py $_.FullName
>>     $output | Out-File -FilePath $outputFile
```

### Problem 5

```shell
Get-ChildItem *.fna | ForEach-Object {
>>     $file = $_.FullName
>>     $outputFile = "$($_.BaseName)_output5.txt"
>>     python.exe Genome5.py -f $file -l 100 | Out-File $outputFile
```

### Problem 6
```shell
Get-ChildItem *.fna | ForEach-Object {
>>     $file = $_.FullName
>>     $outputFile = "$($_.BaseName)_output6.txt"
>>     python.exe Genome6.py $file 20 | Out-File $outputFile
```
