---
title: "Replication_Trevor_Harrington"
author: "Trevor Harrington"
format: html
editor: visual
execute: 
  keep-md: TRUE
  code-fold: TRUE
---





# Week 1 -- From Statistical Analysis to Bioinformatics

## Sampling, Processing, and Analyzing Genomic data

The goal of this notebook is to complete a series of challenges that will lay the foundation for using RStudio to analyze genomic data. By learning how to effectively store objects in variables, build functions to analyze genetic data, and implement "for" loops when necessary, we aim to enhance our ability to analyze and interpret genomic data more efficiently. Through analyzing genomic data with the tools and concepts covered in the course, we will gain a deeper understanding of genetic information and its crucial role in biological processes. Ultimately, our aim is to leave this course with a well-documented online environment that showcases our accomplishments in analyzing genomic data

------------------------------------------------------------------------

### Challenge 1: Generate a Random List of DNA Nucleotides


::: {.cell}

```{.r .cell-code}
# create a new dataframe with two columns
nt_abr <- c("A", "T", "G", "C")

# generate a new value representing the number of nucleotides being randomly created for this dataset
nt_num <- 1500

set.seed(215) # for reproducibility, using the same seed each time means the "randomness" can be recreated.
nt_sample <- sample(nt_abr, size = nt_num, replace = TRUE)

paste(nt_sample, collapse = "")
```

::: {.cell-output .cell-output-stdout}
```
[1] "CGGAACTCCCAACGCCTTGATTCCCGAGTTCTAAGCCGGATCATTGTGGTTTTTGATTGAGAGTCAATCTCAAACGACGTAAGTAGTGTGCGTTGAGCTCTCGCGGATAGGACTATACCGGACGCGAGTTAAGACTCTGAGACGAAAAATAAGCAGGCCTCTCACTGTCGGTCTTAACTACCCCCACTTCCCGTCGTACATCCACGGTTTCTTAATTCCGTGAACCGTGGTACCATGCCTCACGTATGAAAGAGGTATGAGACGCACACTATCTCCTAGTCACTCGATATAGGCAGGTACCGGACGCTAGTGCATGATTGGGGCAGGTAATTTTCCTCGCCGTTTTGTCGGTTGGACGTTAAGGCGCCCCATACTGGCACCCCAAAGACGCGAATCCAAGCAATTCATCCAGTCGATGGGAAGGCGTATACGCACGCGCGATCTCATATTAAGACGTATCCAGGTACTAACAAAAGCCGTTGGCCCGCAACTTCCGATAGCAACTGAGACAACCATCCGTGGGTATACAATTCATCTGGCCTGCTTTTTCCAGTGCATGGGGGGGCGGGGAATTACTAGTGCCTGCACGCCAGCTTCTTGGTACCCCCGGCTATGGTTTCTTACGAAACATACATCAGTGCGTTACCTTGGCTAAACGGTTTCGTGAAGGACGTGGCGTGGCAAGTGCGCGGGTTACATCCGGGTAGACCGAGCCACAAGAAATAGGTAAACCGGGAATCAGCGTGGTATACGCATAGTCTGTATCCTTCGGGGGGTGATCATTGAACCTTCAGCACGCTCGAGCAGTGCATTGGGTTAGTCCCGGGCTTCTCCATTCCGCAACGGGACTGGTTCCACGGGACGTTACGATGATATCGTGGGAGGCCAACAAGCGCTTGTGACAAAGTTCTGCGGCTGAAACTCGCCATGCGTCTCCTCGCTACCCGTTACACCGGCAAAGCCTAAGAGTTATTTCACTCACCGATCTTCAGACTCTTTAAACGCACTTCTTAATAGCATCCTATTCGGACGGCAATTTTATTCTAACTATTATCCAGGCTCTTAGCATCCCCGAGTTGTTCATGTGTCTCTGAAACAAGTCAAGTCAGTCATGGCTTGCGCCAGGTGAAAAACGAACTTTCCCTTAGTGTCATTAAGCGCAGCTGGACGTGGCCGATCCACGTTCGTTTGGCGAGTAGAACCCAAACTGTCGCATTAACTAGTATCATTATAAGGTGAATGGGGTAGTTATGATCGGTGTTATAAGTACGGAAGGAGCGCTGGGTGAAACTAGGGTTGTACCAAGTTGCCCCACTATGGGCAAAGAGGCTATACGCGTATAGAGTTAGACAGCTGGTAGCGAAGGGTATAATGTCCCCCTGGTAGGAAGATTTCCGTATTCCCACTATATTGTGGGCTACTATTGGTTAAAGTAGCTGGTGGCTAATGCCAAAACCACCGTACGCTCAGGCAAAGTAGGATTAGAGAGATGCTGAAACT"
```
:::
:::


------------------------------------------------------------------------

### Challenge 2: Generate random nucleotides sequences 15 characters long


::: {.cell}

```{.r .cell-code}
# create a new dataframe with two columns
nt_abr <- c("A", "C", "G", "T")

set.seed(777) # for reproducibility
rnd_genome <- sample(nt_sample, size = 15, replace = TRUE)

paste(rnd_genome, collapse = "")
```

::: {.cell-output .cell-output-stdout}
```
[1] "CGCGGAGACATGGTG"
```
:::
:::


------------------------------------------------------------------------

### Challenge 3: generating a randomized 1500 nucleotide-long dataset


::: {.cell}

```{.r .cell-code}
# create a new dataframe with two columns
nt_abr <- c("A", "T", "G", "C")

# generate a new value representing the number of nucleotides being randomly created for this dataset
nt_num <- 1500

set.seed(215) # for reproducibility, using the same seed each time means the "randomness" can be recreated.
nt_sample <- sample(nt_abr, size = nt_num, replace = TRUE)

paste(nt_sample, collapse = "")
```

::: {.cell-output .cell-output-stdout}
```
[1] "CGGAACTCCCAACGCCTTGATTCCCGAGTTCTAAGCCGGATCATTGTGGTTTTTGATTGAGAGTCAATCTCAAACGACGTAAGTAGTGTGCGTTGAGCTCTCGCGGATAGGACTATACCGGACGCGAGTTAAGACTCTGAGACGAAAAATAAGCAGGCCTCTCACTGTCGGTCTTAACTACCCCCACTTCCCGTCGTACATCCACGGTTTCTTAATTCCGTGAACCGTGGTACCATGCCTCACGTATGAAAGAGGTATGAGACGCACACTATCTCCTAGTCACTCGATATAGGCAGGTACCGGACGCTAGTGCATGATTGGGGCAGGTAATTTTCCTCGCCGTTTTGTCGGTTGGACGTTAAGGCGCCCCATACTGGCACCCCAAAGACGCGAATCCAAGCAATTCATCCAGTCGATGGGAAGGCGTATACGCACGCGCGATCTCATATTAAGACGTATCCAGGTACTAACAAAAGCCGTTGGCCCGCAACTTCCGATAGCAACTGAGACAACCATCCGTGGGTATACAATTCATCTGGCCTGCTTTTTCCAGTGCATGGGGGGGCGGGGAATTACTAGTGCCTGCACGCCAGCTTCTTGGTACCCCCGGCTATGGTTTCTTACGAAACATACATCAGTGCGTTACCTTGGCTAAACGGTTTCGTGAAGGACGTGGCGTGGCAAGTGCGCGGGTTACATCCGGGTAGACCGAGCCACAAGAAATAGGTAAACCGGGAATCAGCGTGGTATACGCATAGTCTGTATCCTTCGGGGGGTGATCATTGAACCTTCAGCACGCTCGAGCAGTGCATTGGGTTAGTCCCGGGCTTCTCCATTCCGCAACGGGACTGGTTCCACGGGACGTTACGATGATATCGTGGGAGGCCAACAAGCGCTTGTGACAAAGTTCTGCGGCTGAAACTCGCCATGCGTCTCCTCGCTACCCGTTACACCGGCAAAGCCTAAGAGTTATTTCACTCACCGATCTTCAGACTCTTTAAACGCACTTCTTAATAGCATCCTATTCGGACGGCAATTTTATTCTAACTATTATCCAGGCTCTTAGCATCCCCGAGTTGTTCATGTGTCTCTGAAACAAGTCAAGTCAGTCATGGCTTGCGCCAGGTGAAAAACGAACTTTCCCTTAGTGTCATTAAGCGCAGCTGGACGTGGCCGATCCACGTTCGTTTGGCGAGTAGAACCCAAACTGTCGCATTAACTAGTATCATTATAAGGTGAATGGGGTAGTTATGATCGGTGTTATAAGTACGGAAGGAGCGCTGGGTGAAACTAGGGTTGTACCAAGTTGCCCCACTATGGGCAAAGAGGCTATACGCGTATAGAGTTAGACAGCTGGTAGCGAAGGGTATAATGTCCCCCTGGTAGGAAGATTTCCGTATTCCCACTATATTGTGGGCTACTATTGGTTAAAGTAGCTGGTGGCTAATGCCAAAACCACCGTACGCTCAGGCAAAGTAGGATTAGAGAGATGCTGAAACT"
```
:::
:::


------------------------------------------------------------------------

### Challenge 4: Counting each DNA Nucleotide in a Randomized Dataset and using Loop()

#### Generating a 100-nucleotide strand and counting each frequency


::: {.cell}

```{.r .cell-code}
# create a new dataframe with two columns
nt_abr <- c("A", "C", "G", "T")

set.seed(777) # for reproducibility
rnd_genome <- sample(nt_sample, size = 100, replace = TRUE)
table(rnd_genome)
```

::: {.cell-output .cell-output-stdout}
```
rnd_genome
 A  C  G  T 
20 33 30 17 
```
:::
:::


-   This distribution heavily favors C va

#### Using Loop() function to repeat an operation


::: {.cell}

```{.r .cell-code}
rnd_sum <- 0

for(i in 1:10){
  rnd_sum <- rnd_sum + i
  print(rnd_sum)
}
```

::: {.cell-output .cell-output-stdout}
```
[1] 1
[1] 3
[1] 6
[1] 10
[1] 15
[1] 21
[1] 28
[1] 36
[1] 45
[1] 55
```
:::
:::


------------------------------------------------------------------------

### Challenge 5: Using loop() functions with while/for staements to generate a


::: {.cell}

```{.r .cell-code}
#initialize a new dataframe with a value of 1
myProduct <- 1

# Use "for loop" to iterate (multiply) myProduct by values of j(1-15)
##kind of confusing, but the iteration variable "j" 
for(j in 1:15) {
  myProduct <- myProduct * j   # Update myProduct by multiplying with j
  print(myProduct)            # Print the updated value of myProduct
}
```

::: {.cell-output .cell-output-stdout}
```
[1] 1
[1] 2
[1] 6
[1] 24
[1] 120
[1] 720
[1] 5040
[1] 40320
[1] 362880
[1] 3628800
[1] 39916800
[1] 479001600
[1] 6227020800
[1] 87178291200
[1] 1.307674e+12
```
:::
:::


------------------------------------------------------------------------

### Challenge 6: Combining sample(), paste() and loop() functions to produce a set of individual nucleotides

#### 

Generating Random Genomic Data in a collapsed string


::: {.cell}

```{.r .cell-code}
nt <- c("A", "C", "G", "T")
genomeLength <- 10

set.seed(215)
rnd_genome <- paste(
  sample(nt, size = genomeLength, replace = TRUE),
                   collapse = "")
print(rnd_genome)
```

::: {.cell-output .cell-output-stdout}
```
[1] "TGGAATCTTT"
```
:::
:::


#### Counting Adenine using loop() function with If/else statements


::: {.cell}

```{.r .cell-code}
# "For loop" is used to perform iterations 1:n where n (nchar) refers to the number of characters (variables, nucleotides) in the rnd_genome dataset
for(i in 1:nchar(rnd_genome)){
  # String subset is checked for each iteration, start = i and end = i means this substring will only extract a single character, i, for each iteration
  ## == "A" checks each substring, and only extracts A for each nucleotide in the datafrom
  if(str_sub(rnd_genome, start = i, end = i) == "A"){
    print(str_sub(rnd_genome, start = i, end = i))
  }
}
```

::: {.cell-output .cell-output-stdout}
```
[1] "A"
[1] "A"
```
:::
:::


#### Modifying Output to sum the Variables Extracted from for Loop() function


::: {.cell}

```{.r .cell-code}
sum_A <- 0  # Create a new variable to store sum of "A" characters extracted from the loop

for(i in 1:nchar(rnd_genome)) {
  if(str_sub(rnd_genome, start = i, end = i) == "A") {
    sum_A <- sum_A + 1   # Adds 1 to sum_A variable when a iterations passes the if "A" statement
  }
}

cat("Total number of Adenisine Nucleotides: ", sum_A)  # Print the final value of sum_A using cat() "concatenate and print" function to merge the output of sum_A
```

::: {.cell-output .cell-output-stdout}
```
Total number of Adenisine Nucleotides:  2
```
:::

```{.r .cell-code}
## Without cat(), print() produces two lines (1) and (2)
```
:::


### Challenge 7: 

# Understanding the biology -- DNA Replication Lecture 

## Useful Terminology

-   Scales of genetic data; genome -\> Chromosome -\> Gene -\> Nucleotide

    -   Humans have 23 pairs of chromosomes

    -   Bacteria have a single looping strand containing all of their genetic data

        -   Bacteria also contain small, circular DNA molecules called **Plasmids**, used for storing additional data and conjugation/cloning

-   **Complementary base pairing:** A & T + C & G

    -   antiparellel strands: DNA always read from **5' -\> 3'**

    -   DNA strand = phosphate group + hydroxyl (sugar) group + base group (A,C,G,T)

        -   Sugar-Phosphate backbone: creates the 5' -\> 3' rule (placement of carbon in sugar)

        -   base pairs held together by matching **Hydrogen bonds.**

-   **DNA Polymerase:** Along with a number of other enzymes used to read and create a RNA copy of the DNA.

    ## **Origins of Replication**

    shorts sequences of genes that can be read to initiate the reading of a gene.

-   Many origins of replication in the human genome, where genes are constantly being replicated

-   Bacteria (tend to) only have a single origin of replication that travels in both directions

    -   **DnaA Boxes:** specific bacterial marker that identifies the origin for enzymes to begin replication

        -   specific sequence used to initiate a DnaA box
