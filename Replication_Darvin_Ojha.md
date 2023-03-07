---
title: "Starter Notebook"
author: "You, Scientist"
format: html
execute:
  keep-md: true
---





## Challenge 1 : Listing four nucleotides


::: {.cell}

```{.r .cell-code}
nucl_names <- c("A", "T", "G", "C")
nucl_names #nucl_names represents the names of the nucleotides. 
```

::: {.cell-output .cell-output-stdout}
```
[1] "A" "T" "G" "C"
```
:::
:::


## Challenge 2 : Creating a random string of 15 nucleotides


::: {.cell}

```{.r .cell-code}
randGenome <- sample(nucl_names, size = 15, replace = TRUE)
paste(randGenome, collapse = "")
```

::: {.cell-output .cell-output-stdout}
```
[1] "GTATGTGCCAATCCC"
```
:::
:::


### Challenge 3: Generate a random genome which is 1500 nucleotides long dataset


::: {.cell}

```{.r .cell-code}
set.seed(215)
genomeLength <- 1500
randGenome <- sample(nucl_names, size = 1500, replace = TRUE)
paste(randGenome, collapse = "")
```

::: {.cell-output .cell-output-stdout}
```
[1] "CGGAACTCCCAACGCCTTGATTCCCGAGTTCTAAGCCGGATCATTGTGGTTTTTGATTGAGAGTCAATCTCAAACGACGTAAGTAGTGTGCGTTGAGCTCTCGCGGATAGGACTATACCGGACGCGAGTTAAGACTCTGAGACGAAAAATAAGCAGGCCTCTCACTGTCGGTCTTAACTACCCCCACTTCCCGTCGTACATCCACGGTTTCTTAATTCCGTGAACCGTGGTACCATGCCTCACGTATGAAAGAGGTATGAGACGCACACTATCTCCTAGTCACTCGATATAGGCAGGTACCGGACGCTAGTGCATGATTGGGGCAGGTAATTTTCCTCGCCGTTTTGTCGGTTGGACGTTAAGGCGCCCCATACTGGCACCCCAAAGACGCGAATCCAAGCAATTCATCCAGTCGATGGGAAGGCGTATACGCACGCGCGATCTCATATTAAGACGTATCCAGGTACTAACAAAAGCCGTTGGCCCGCAACTTCCGATAGCAACTGAGACAACCATCCGTGGGTATACAATTCATCTGGCCTGCTTTTTCCAGTGCATGGGGGGGCGGGGAATTACTAGTGCCTGCACGCCAGCTTCTTGGTACCCCCGGCTATGGTTTCTTACGAAACATACATCAGTGCGTTACCTTGGCTAAACGGTTTCGTGAAGGACGTGGCGTGGCAAGTGCGCGGGTTACATCCGGGTAGACCGAGCCACAAGAAATAGGTAAACCGGGAATCAGCGTGGTATACGCATAGTCTGTATCCTTCGGGGGGTGATCATTGAACCTTCAGCACGCTCGAGCAGTGCATTGGGTTAGTCCCGGGCTTCTCCATTCCGCAACGGGACTGGTTCCACGGGACGTTACGATGATATCGTGGGAGGCCAACAAGCGCTTGTGACAAAGTTCTGCGGCTGAAACTCGCCATGCGTCTCCTCGCTACCCGTTACACCGGCAAAGCCTAAGAGTTATTTCACTCACCGATCTTCAGACTCTTTAAACGCACTTCTTAATAGCATCCTATTCGGACGGCAATTTTATTCTAACTATTATCCAGGCTCTTAGCATCCCCGAGTTGTTCATGTGTCTCTGAAACAAGTCAAGTCAGTCATGGCTTGCGCCAGGTGAAAAACGAACTTTCCCTTAGTGTCATTAAGCGCAGCTGGACGTGGCCGATCCACGTTCGTTTGGCGAGTAGAACCCAAACTGTCGCATTAACTAGTATCATTATAAGGTGAATGGGGTAGTTATGATCGGTGTTATAAGTACGGAAGGAGCGCTGGGTGAAACTAGGGTTGTACCAAGTTGCCCCACTATGGGCAAAGAGGCTATACGCGTATAGAGTTAGACAGCTGGTAGCGAAGGGTATAATGTCCCCCTGGTAGGAAGATTTCCGTATTCCCACTATATTGTGGGCTACTATTGGTTAAAGTAGCTGGTGGCTAATGCCAAAACCACCGTACGCTCAGGCAAAGTAGGATTAGAGAGATGCTGAAACT"
```
:::
:::


### Challenge 4: Generating random genome consisting of 100 nucleotides


::: {.cell}

```{.r .cell-code}
set.seed(215)
genomeLength <- 100
randGenome <- sample(nucl_names, size = 100, replace = TRUE)
table(randGenome)
```

::: {.cell-output .cell-output-stdout}
```
randGenome
 A  C  G  T 
23 23 25 29 
```
:::
:::


### Challenge 5: Writing my own loop 


::: {.cell}

```{.r .cell-code}
mySum <- 0

for(i in 1:10){
  mySum <- mySum + i
  print(mySum)
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


### Challenge 6 : 

### Generating a random genome substring consisting of 10 nucleotides using sample(),paste() and loop


::: {.cell}

```{.r .cell-code}
nucleotides <- c("A", "C", "G", "T")
genomeLength <- 10

randGenome <- paste(
  sample(nucleotides, size = genomeLength, replace = TRUE),
                   collapse = "")
print(randGenome)
```

::: {.cell-output .cell-output-stdout}
```
[1] "CTGTGGACAG"
```
:::
:::


### Using for loop to print each nucleotides in my random sample


::: {.cell}

```{.r .cell-code}
for(i in 1:nchar(randGenome)){
  
    print(str_sub(randGenome, start = i, end = i))
  }
```

::: {.cell-output .cell-output-stdout}
```
[1] "C"
[1] "T"
[1] "G"
[1] "T"
[1] "G"
[1] "G"
[1] "A"
[1] "C"
[1] "A"
[1] "G"
```
:::
:::


### Challenge 7 :  Counting the number of occurrences of Adenine (A) in randGenome.


::: {.cell}

```{.r .cell-code}
sum_A <- 0

for(i in 1:nchar(randGenome)){
  if(str_sub(randGenome, start = i, end = i) == "A")
    sum_A <- sum_A + 1
    {
    print(sum_A)
  }
}
```

::: {.cell-output .cell-output-stdout}
```
[1] 0
[1] 0
[1] 0
[1] 0
[1] 0
[1] 0
[1] 1
[1] 1
[1] 2
[1] 2
```
:::
:::


### Challenge 8 : Counting the frequencies of each of the four individual nucleotides.


::: {.cell}

:::
