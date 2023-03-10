---
title: "Replication_Trevor_Harrington"
author: "Trevor Harrington"
format: html
editor: visual
execute: 
  keep-md: TRUE
  code-fold: TRUE
---





# Chapter 1 -- From Statistical Analysis to Bioinformatics

## Introduction: Sampling, Processing, and Analyzing Genomic data

The goal of this notebook is to complete a series of challenges that will lay the foundation for using RStudio to analyze genomic data. By learning how to effectively store objects in variables, build functions to analyze genetic data, and implement "for" loops when necessary, we aim to enhance our ability to analyze and interpret genomic data more efficiently. Through analyzing genomic data with the tools and concepts covered in the course, we will gain a deeper understanding of genetic information and its crucial role in biological processes. Ultimately, our aim is to leave this course with a well-documented online environment that showcases our accomplishments in analyzing genomic data

## Objectives:

-   Use the arrow (\<-) operator to store objects in variables.

-   Build functions to process and analyze genetic data.

-   Create a list of nucleotides

-   Use set.seed() to set a seed for a random number generator to make your results reproducible over time and across different machines (and different researchers too).

-   Use R's sample() function to create a random "genome" for the purpose of testing functions.

-   Use paste() with the parameter setting collapse = "" to collapse your random genome into a single long string of nucleotides.

-   Learn about `for loops` and implement them when necessary.

-   Read in a real bacterial genome from a text file.

-   Use programmatic flow control in the form of if/else if/else statements to run code only when a particular condition is satisfied.

-   Solve a bioinformatics problem on the Rosalind platform.

------------------------------------------------------------------------

### Challenge 1: Generate a Random List of DNA Nucleotides


::: {.cell}

```{.r .cell-code}
# create a new variable containing a vector of four characters for the four nucleotide abbreviations (Adenine, Cytosine, Guanine, Thymine) 
## Vector is  a one-dimensional array-like object that can contain values of the same basic type of data (numerical,character,logical)
nt_names <- c("A", "T", "G", "C")

nt_names # simply displays the vectors without organizing
```

::: {.cell-output .cell-output-stdout}
```
[1] "A" "T" "G" "C"
```
:::
:::


------------------------------------------------------------------------

### Challenge 2: Generate 15 Character String of Nucleotides


::: {.cell}

```{.r .cell-code}
set.seed(215)

# `sample` function used to generate a random sample of nucleotides from the nucleotides (nt_name) vector. The `size` argument specifies the number of samples to generate (15)  
##`replace = TRUE` indicates sample should be replaced in the pool after being selected (same nucleotide can appear multiple times in the sample)
## nt_sample calls on vector generated in challenge 1 and stored in the environment
rnd_genome <- sample(nt_names, size = 15, replace = TRUE)

paste(rnd_genome, collapse = "")
```

::: {.cell-output .cell-output-stdout}
```
[1] "CGGAACTCCCAACGC"
```
:::
:::


------------------------------------------------------------------------

### Challenge 3: Generating a Randomized 1500 Nucleotide-Long Dataset


::: {.cell}

```{.r .cell-code}
# create a new variable containing a vector of four characters for the four nucleotide abbreviations (Adenine, Cytosine, Guanine, Thymine) 
## Vector is  a one-dimensional array-like object that can contain values of the same basic type of data (numerical,character,logical)
nt_names <- c("A", "T", "G", "C")

# generate a new value representing the number of nucleotides that will be randomly generated for this dataset
##This step is not required (numerical value can be used in 'sample') but it is better to call on a value then hard code one into a function. More flexible when df is larger
nt_num <- 1500


# for reproducibility, using the same seed each time means the "randomness" can be recreated.
set.seed(215)


#'sample' function randomly picks variables from nt_names vector to generate a new vector n length, in this case n = nt_num (1500)
nt_sample <- sample(nt_names, size = nt_num, replace = TRUE)

# 'paste' combines the individual nucleotides into a single string rather then a list of individual letters
# 'collapse' removes the seperator generated between each nucleotide as a result of the 'sample' function interavting with the nucleotides.
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
# Set a seed value to ensure that the random numbers generated are reproducible
set.seed(215)

# Generate a random genome sequence of length 100 by sampling (with replacement) from the nucleotide abbreviation vector still in the environment
rnd_genome <- sample(nt_names, size = 100, replace = TRUE)

# 'table' function counts the frequency of each vector (nucleotide) in the random genome sequence and generates a table of the results 
table(rnd_genome)
```

::: {.cell-output .cell-output-stdout}
```
rnd_genome
 A  C  G  T 
23 23 25 29 
```
:::
:::


#### Using Loop() function to repeat an operation

-   `for loop` functions are widely used and incredibly versatile for repeating an operation with minimal inputs.


::: {.cell}

```{.r .cell-code}
# Create 'loop_result' variable used to store the values being added in each iteration of the loop.
## loop is initialized to 0 (or whatever value entered for 'loop_result') and then incremented by the current value of i (and the function applied to the loop) in each iteration of the loop
loop_result <- 0

# for loop is a control flow statement in programming that allows you to execute a block of code repeatedly.
## The loop consists of an initialization statement, a test condition, an update statement, and a loop body. In this code
## Here, a loop is used to iterate over a sequence of values from 1 to 10, adding the value 'i' for each iteration.
for(i in 1:10){
  loop_result <- loop_result + i
  print(loop_result)
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

### Challenge 5: Using loop() functions


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

#### Generating Random Genomic Data in a collapsed string


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
for(i in 1:nchar(rnd_genome))
  {
    print(str_sub(rnd_genome, start = i, end = i))
  }
```

::: {.cell-output .cell-output-stdout}
```
[1] "T"
[1] "G"
[1] "G"
[1] "A"
[1] "A"
[1] "T"
[1] "C"
[1] "T"
[1] "T"
[1] "T"
```
:::
:::


### Challenge 7: Modifying Output to sum the For Loop() Results Output


::: {.cell}

```{.r .cell-code}
# Set a variable called `sum_A` equal to 0
sum_A <- 0

# Use a for loop to iterate over each character in the `rnd_genome` variable
for (i in 1:nchar(rnd_genome)) {
  # # String subset is checked for each iteration, start = i and end = i means this substring will only extract a single character, i, for each iteration
  ## == "A" checks each substring, and only extracts A for each nucleotide in the datafrom
  if (str_sub(rnd_genome, start = i, end = i) == "A") {
    # If it is, add 1 to the `sum_A` variable
    sum_A <- sum_A + 1
  }
}

# Create a named vector called `result` with one element, the number of adenine bases found
result <- c("Adenine:" = sum_A)

# Use the `print` function to output the `result` vector to the console
print(result)
```

::: {.cell-output .cell-output-stdout}
```
Adenine: 
       2 
```
:::
:::

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
## Without cat(), print() produces two lines
```
:::


### Challenge 8: Adding the remaining nucleotides to the count


::: {.cell}

```{.r .cell-code}
# Initialize count variables for each nucleotide to 0
sum_A <- 0
sum_C <- 0
sum_G <- 0
sum_T <- 0

# Loop through each character in the random genome sequence
for(i in 1:nchar(rnd_genome)) {
  # Check if the character is "A", "C", "G", or "T", and increment the corresponding count variable if it is
  ##'str_sub' substring function allows each vector in the variable to be checked and summed independent of the other vectors in the string
  ##'start = i' and 'end = i' is starting at position i and ending at position i for each iteration extracting the "i-th" character of the string, which is then compared to the nucleotide letters A, C, G, and T and added to the count.
  if(str_sub(rnd_genome, start = i, end = i) == "A") {
    sum_A <- sum_A + 1}
  if(str_sub(rnd_genome, start = i, end = i) == "C") {
    sum_C <- sum_C + 1}
  if(str_sub(rnd_genome, start = i, end = i) == "G") {
    sum_G <- sum_G + 1}
  if(str_sub(rnd_genome, start = i, end = i) == "T") {
    sum_T <- sum_T + 1}
}

# Create a named vector of the counts for each nucleotide
result <- c("Adenine:" = sum_A, "Cytosine:" = sum_C, "Guanine:" = sum_G, "Thymine:" = sum_T)

# Print the results
print(result)
```

::: {.cell-output .cell-output-stdout}
```
 Adenine: Cytosine:  Guanine:  Thymine: 
        2         1         2         5 
```
:::
:::


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

------------------------------------------------------------------------


::: {.cell}

```{.r .cell-code}
vib_c <- scan("C:/School/23SPDAY/Biostatistics/Group_Project/Genomes/VibrioCholerae.txt", what = "character", sep = NULL)
```
:::

::: {.cell}

```{.r .cell-code}
sum_A <- 0
sum_C <- 0
sum_G <- 0
sum_T <- 0

for(i in 1:nchar(vib_c)) {
  if(str_sub(vib_c, start = i, end = i) == "A") {sum_A <- sum_A + 1}
  if(str_sub(vib_c, start = i, end = i) == "C") {sum_C <- sum_C + 1}
  if(str_sub(vib_c, start = i, end = i) == "G") {sum_G <- sum_G + 1}
  if(str_sub(vib_c, start = i, end = i) == "T") {sum_T <- sum_T + 1}
}

result <- c("Adenine:" = sum_A, "Cytosine:" = sum_C, "Guanine:" = sum_G, "Thymine:" = sum_T)
print(result)
```

::: {.cell-output .cell-output-stdout}
```
 Adenine: Cytosine:  Guanine:  Thymine: 
   293942    263573    256024    294711 
```
:::
:::


------------------------------------------------------------------------

### Rosalind Testing Data


::: {.cell}

```{.r .cell-code}
rosalind_string <- "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
```
:::

::: {.cell}

```{.r .cell-code}
sum_A <- 0
sum_C <- 0
sum_G <- 0
sum_T <- 0

for(i in 1:nchar(rosalind_string)) {
  if(str_sub(rosalind_string, start = i, end = i) == "A") {sum_A <- sum_A + 1}
  if(str_sub(rosalind_string, start = i, end = i) == "C") {sum_C <- sum_C + 1}
  if(str_sub(rosalind_string, start = i, end = i) == "G") {sum_G <- sum_G + 1}
  if(str_sub(rosalind_string, start = i, end = i) == "T") {sum_T <- sum_T + 1}
}

paste(sum_A,sum_C,sum_G,sum_T)
```

::: {.cell-output .cell-output-stdout}
```
[1] "20 12 17 21"
```
:::
:::


-   This result matches the key provided on the website, so I assume that this code is effective!

------------------------------------------------------------------------

# Chapter 2: Bioinformatic Tools pt. 2

## Introduction: Improving Tools and Efficiency

-   Genomes are large, and it is important that the tools for analyzing them are used properly to both validate that the results are accurate, and that as few unnecessary steps as possible are made since they will be time consuming when working with a large string of nucleotides.

## Objectives:

-   Learn to eliminate "single-use code" by encapsulating code in functions which can be easily applied to different genomes.

-   Expand your ability to count the frequency of each individual nucleotide within a genome to an ability to count the appearance of patterns in the genome.

-   Build and execute functions to find "hidden messages" within the genome by identifying patterns appearing much more often than they would be expected to if a genome included nucleotides at random.

-   Exploit transcription errors to narrow a search region for the replication origin.

------------------------------------------------------------------------

### Challenge 0: Generating Functions for Functions


::: {.cell}

```{.r .cell-code}
# Functions allow code to be ran without having to enter or copy and paste the code into a new block
## within function(), the variables that can be altered are added (and pop up when entering in the function)
## nucleotide = "A" means IF no nucleotide is supplied when using the function, it defaults to "A"
nucleotide_frequency <- function(genomeString, nucleotide = "A"){
  count <- 0 # generate a new variable for the function to write to
  
  ## for each iteration of i from 1 to the end of 'genomestring', if character in the string is the specificied nucleotide, count is updated by +1
## {} are layers to the formula, which are only run in each iteration if the previous step is satisfied.
  for(i in 1:nchar(genomeString)){
    if(str_sub(genomeString, start = i, end = i) == nucleotide){
      count <- count + 1
    }
  }
  return(count)
}

nucleotide_frequency("ACTTGCGGGTATCGAG", "G")
```

::: {.cell-output .cell-output-stdout}
```
[1] 6
```
:::
:::


------------------------------------------------------------------------

### Challenge 1: Using the Function Function on a Randomized Genome


::: {.cell}

```{.r .cell-code}
nt_sample <- sample(nt_names, size = 2000, replace = TRUE)
nt_sample <- paste(nt_sample, collapse = "")
```
:::

::: {.cell}

```{.r .cell-code}
nucleotide_frequency(nt_sample, "C")
```

::: {.cell-output .cell-output-stdout}
```
[1] 496
```
:::
:::


-   23.8% of the nucleotides in this randomized dataset are cytosines

    -   code is random, but this character

### Generating a function to create randomized genome


::: {.cell}

```{.r .cell-code}
nt_names <- c("A", "C", "G", "T")

rnd_genome <- function(k){
  set.seed(215)
  rnd_string <- ""
  
  for (i in 1:k){
    rnd_string <- paste(rnd_string, sample(nt_names,1, replace = TRUE), sep="")
  }
  return(rnd_string)
}
```
:::

::: {.cell}

```{.r .cell-code}
rnd_genome(9)
```

::: {.cell-output .cell-output-stdout}
```
[1] "TGGAATCTT"
```
:::
:::


#### From single-vector multiple-vector Substrings

-   `str_sub` or substrings, allow the function to extract characters based on a starting and and ending point, which can be valuable for extracting vectors of a specific length, or with specific values. This function will be valuable for genomic data to identify patterns and extract vectors with specific lengths.


::: {.cell}

```{.r .cell-code}
string_sample <- c()
# subset strings extract a value based on the 
string_sample <- string_sample %>% 
  append(str_sub(nt_sample, start = 1000, end = 1005))

# Negative vector position indicated the position to start and end from is being measured starting at the end of the dataset
# When data from variables is added to a substring, they are removed and cannot be added to a new substring that overlaps it
str_sub(nt_sample, start = -1000, end = -995)
```

::: {.cell-output .cell-output-stdout}
```
[1] "TTAATA"
```
:::
:::


#### Appending a new list containing a substring of nucleotides


::: {.cell}

```{.r .cell-code}
# 'append' is useful for generating a reusable variable containing a vector that can be called back for later tests
string_sample <- c()
# subset strings exctract a value based on the 
string_sample <- string_sample %>% 
  append(str_sub(nt_sample, start = 1000, end = 1005))

string_sample
```

::: {.cell-output .cell-output-stdout}
```
[1] "CTTAAT"
```
:::
:::


#### Using multiple-vector substrings to generate paired nucleotides in a dataset


::: {.cell}

```{.r .cell-code}
generate_2_mer <- function(string_sample) {
  list_codon <- c()

  
  # this function will run from 1 to the end of the supplied genome sting (-1 prevents the function from adding a variable that does not contain 2 variables)
  for(i in 1:(nchar(string_sample) - 1)){
  list_codon <- list_codon %>%
  append(str_sub(string_sample, start = i, end = i + 1))
    }
  return(list_codon)
}

generate_2_mer(string_sample)
```

::: {.cell-output .cell-output-stdout}
```
[1] "CT" "TT" "TA" "AA" "AT"
```
:::
:::


### Using multiple-vector substrings to generate codons in a dataset


::: {.cell}

```{.r .cell-code}
generate_codons <- function(string_sample){
  list_codon <- c()

# for each iteration i in the sequance 1 through all charictars in the provided string 
  ## -2 to prevent a partial codon
  ## by = 3 to make the function shift 3 rather then shifting 1 and counting 3 for each iteration
  
  for (i in seq(1, nchar(string_sample) - 2, by = 3)) {
    list_codon <- list_codon %>%
      ##append adds an additional vector to the variable, without clearing what was generated previously. end = i + 2 means the function will count 3 out from the current variable being iterated, which jumps 3 each iteration. 
      
      append(str_sub(string_sample, start = i, end = i + 2))
  }
  
  return(list_codon)
}

generate_codons(rnd_genome(200))
```

::: {.cell-output .cell-output-stdout}
```
 [1] "TGG" "AAT" "CTT" "TAA" "TGT" "TCC" "GAC" "CTT" "TGA" "GCC" "TCA" "AGT"
[13] "TGG" "ACT" "ACC" "GCG" "GCC" "CCC" "GAC" "CGA" "GAG" "CTA" "ACT" "CTA"
[25] "AAT" "GAT" "GCA" "AGC" "AGC" "GCG" "TGC" "CGA" "GTC" "TCT" "GTG" "GAC"
[37] "AGG" "ATC" "ACA" "TTG" "GAT" "GTG" "AGC" "CAA" "GAT" "CTC" "GAG" "ATG"
[49] "AAA" "AAC" "AAG" "TAG" "GTT" "CTC" "TAT" "CGC" "TGG" "CTC" "CAA" "TCA"
[61] "TTT" "TTA" "TCC" "TTT" "GCT" "GCA"
```
:::
:::

::: {.cell}

```{.r .cell-code}
generate_k_mer <- function(string, k = 3) {
  list_codon <- c()

  for (i in seq(1, nchar(string) - k + 1, by = k)) {
    list_codon <- list_codon %>%
      append(str_sub(string, start = i, end = i + k - 1))
  }
  
  return(list_codon)
}

generate_k_mer(rnd_genome(9))
```

::: {.cell-output .cell-output-stdout}
```
[1] "TGG" "AAT" "CTT"
```
:::
:::


-   This function can

### Using functions to generate genome, extract strings of length k, and find the frequency of nucleotides 


::: {.cell}

```{.r .cell-code}
generate_k_mer(rnd_genome(1000),6)
```

::: {.cell-output .cell-output-stdout}
```
  [1] "TGGAAT" "CTTTAA" "TGTTCC" "GACCTT" "TGAGCC" "TCAAGT" "TGGACT" "ACCGCG"
  [9] "GCCCCC" "GACCGA" "GAGCTA" "ACTCTA" "AATGAT" "GCAAGC" "AGCGCG" "TGCCGA"
 [17] "GTCTCT" "GTGGAC" "AGGATC" "ACATTG" "GATGTG" "AGCCAA" "GATCTC" "GAGATG"
 [25] "AAAAAC" "AAGTAG" "GTTCTC" "TATCGC" "TGGCTC" "CAATCA" "TTTTTA" "TCCTTT"
 [33] "GCTGCA" "TACTTA" "TGGCCC" "TCCAAC" "CTTGCG" "AATTGC" "GGCATT" "ACGTTC"
 [41] "TATGCA" "CGAAAG" "AGGCAC" "GAGATG" "TATATC" "ACTCTT" "CAGCTA" "TCTGAC"
 [49] "ACAGGT" "AGGCAT" "TGGATG" "TCAGCG" "TACGAC" "CGGGGT" "AGGCAA" "CCCCTT"
 [57] "CTGTTG" "CCCCGC" "TGGCCG" "GATGCC" "AAGGTG" "TTTTAC" "ATCGGT" "ATTTTA"
 [65] "AAGATG" "TGAACT" "TAAGTA" "ACCTAC" "TTAGCT" "GACGGG" "AAGGTG" "CACATG"
 [73] "TATGTG" "TGACTC" "TACACC" "AAGATG" "CACTTA" "GGCATC" "AATAAA" "AGTTGC"
 [81] "CGGTTT" "GTAATC" "CTTGAC" "AGTAAT" "CGAGAT" "AATTAC" "TTGCGG" "GCACAT"
 [89] "AACCTA" "CTCGGT" "TCGTCC" "CCCTTA" "GCGTAC" "GGGGGG" "GTGGGG" "AACCAT"
 [97] "CAGCGT" "TCGTAT" "GTTAGT" "CCTCCG" "GCATTT" "TTGGTC" "ACGGCC" "CTCCAT"
[105] "GAAATA" "CATACT" "AGCGTG" "CCATTC" "CGGTCA" "AATGGC" "CCTGCG" "AAGGAT"
[113] "GCGGTG" "CGGTAA" "GCGTGT" "GGGCCA" "TACTTG" "GGCAGA" "TTGAGT" "TATAAG"
[121] "AAACAG" "GCAAAT" "TGGGAA" "CTAGTG" "CGGCAC" "ATGTAC" "AGCTCG" "CACTTC"
[129] "CTGGGG" "GGCGAC" "TACCGA" "ATTCCT" "AGTATG" "TCTGAG" "TAGCGT" "ACCGGG"
[137] "CCAGCT" "TTGGGT" "CCTCTT" "ACCTTG" "TAATGG" "GATCGG" "CCTTAT" "GGGATG"
[145] "CCATGA" "CGACAC" "TGCGGG" "AGGTTA" "ATAAGT" "GTCCGC" "GATAAA" "GCCTCG"
[153] "TGGTCG" "AAATCT" "GTTACG" "TGCTCT" "TCTGTC" "ATTTGC" "CATATT" "GGTAAA"
[161] "GTTCAA" "GAGCCA" "CCCTAT" "CTATTG" "ACTCCT" "AGATCT"
```
:::

```{.r .cell-code}
nucleotide_frequency(rnd_genome(1000),"G")
```

::: {.cell-output .cell-output-stdout}
```
[1] 257
```
:::
:::


-   putting together everything completed in this chapter, the written code can now perform genomic analysis with pre-written functions that are proven to function as intended, and can be utilized for any genome.

### Generating a function to find a specific set of nucleotides


::: {.cell}

```{.r .cell-code}
nt_patterns <- function(string, pattern) {
  nt_matches <- 0
  
  for (i in seq(1,nchar(string))){
    if(str_sub(string, i, i + str_length(pattern)-1) == pattern){
      nt_matches = nt_matches + 1
  }
  }
  return(nt_matches)
}

nt_patterns(rnd_genome(20000), "GACCTT")
```

::: {.cell-output .cell-output-stdout}
```
[1] 9
```
:::
:::

::: {.cell}

```{.r .cell-code}
rosalind_string <- "TTAGTCCCCAGTCCCCAGTCCCCTCCAGTCCCCGCAGTCCCCCAGTCCCCAGTCCCCCAGTCCCCAGTCCCCAGTCCCCAGCCAGTCCCCACCGGTGTGGTAGTCCCCCAAGTCCCCAAGTCCCCAGTGATAACAGTCCCCTTCTCTAAGTCCCCAGTCCCCGAGTCCCCAGTTGAGTCCCCCTAGTCCCCGCCTATAGTCCCCCCACGAGTCCCCTGAAGTCCCCTGAAAGTCCCCCTGACGCAAGTCCCCTAGTCCCCCAGTCCCCAGAAGTCCCCCAGTCCCCTAAGTCCCCTAGAGTCCCCAGTCCCCGAGAGTCCCCTGTAAGTCCCCCTCAGTCCCCGGCTCGAGTCCCCGATGAGTCCCCGAGTCCCCCCGAGTCCCCGGTTAGTCCCCAAGTCCCCGAGTCCCCGAGTCCCCTGAAGTCCCCGAGTCCCCTCGAGTCCCCAGTCCCCGCTAGTCCCCCTTAGTCCCCAGTCCCCGAGTCCCCAGTCCCCAGTCCCCAAGTCCCCCGTGGAGAAGTCCCCGCAGTCCCCAGTCCCCTCGATTAGTCCCCATGCGATAGTCCCCCAGTCCCCTGAGTCCCCAGTCCCCAAGTCCCCGTTAGTCCCCGAGTCCCCAAAATTAGTCCCCGAAGTCCCCCCGTAGTCCCCTGTGAGTCCCCGAGTCCCCAGTCCCCTTACGAGTCCCCGTCCAGTCCCCTGATTATATGAGAGTCCCCTTGGAGTCCCCTAGTCCCCTAAGTCCCCAGAGTGATTCTTAGTCCCCAAGTCCCCAGTCCCCGCTAGTCCCCATAGTCCCCAGTCCCCCGCAGTCCCCCTACCTCAGTCCCCAAGTCCCCCAGTCCCCCAGTCCCCGGTATTAGTCCCCGAAGTCCCCGAGTCCCCATACTCAAGTCCCCCAGTCCCCGGATGGTAGAGTCCCCAGTCCCCTAGTCCCCCCGAGTCCCCAGTCCCCCAGTCCCCACGGCGGCTAAAGTCCCCTAGTCCCCAGTCCCCGATGCAGTCCCCAAGTCCCCAGTCCCC"
```
:::

::: {.cell}

```{.r .cell-code}
nt_patterns(rosalind_string, "AGTCCCCAG")
```

::: {.cell-output .cell-output-stdout}
```
[1] 25
```
:::
:::


-   The code suggests this pattern appeared 25 times, which was the correct answer, meaning the function worked as intended!

------------------------------------------------------------------------

# Chapter 3: Bioinformatic Tools pt. 3
