# protein-annotator
- a genome annotator based on the protein annnotations functions to analyze and annotate your genome using the protein annotations.
- it provides all the analysis done on the basis of the protein anntoations to the genome.
- provide the reference proteins and the assembled genome and it will provide.
  
   - generatemRNA
     ```
     generatemRNA("/home/gaurav/Desktop/final_code_push/multi.gff", 
                        "/home/gaurav/Desktop/final_code_push/multi.fasta", 
                               "/home/gaurav/Desktop/final_code_push/multiout.fasta")
     ```
   - generateCDS
     ```
     generateCDS("/home/gaurav/Desktop/final_code_push/multi.gff", 
                        "/home/gaurav/Desktop/final_code_push/multi.fasta", 
                               "/home/gaurav/Desktop/final_code_push/multiout.fasta")
     ```
   - plotCDS
     ```
     plotCDS("/home/gaurav/Desktop/final_code_push/multi.gff",
                         "/home/gaurav/Desktop/final_code_push/multi.txt")
     ```
     ![codingplotter](https://github.com/sablokgaurav/codingplotter/blob/main/save.png)
   - plotmRNA
     ```
     plotmRNAs("/home/gaurav/Desktop/final_code_push/multi.gff",
                            "/home/gaurav/Desktop/final_code_push/multi.txt")
     ```
   - extractintergenic
     ```
     generateintergenic("/home/gaurav/Desktop/final_code_push/multi.gff",
               "/home/gaurav/Desktop/final_code_push/multi.fasta",
                             "/home/gaurav/Desktop/final_code_push/final.fasta")
     ```
   
## Installation
```bash
$ pip install protein_annotator
```

## Contributing
Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License
`protein_annotator` was created by Gaurav Sablok. It is licensed under the terms of the MIT license. \
Gaurav Sablok \
Academic Staff Member \
Bioinformatics \
Institute for Biochemistry and Biology \
University of Potsdam \
Potsdam,Germany

## Credits

`protein_annotator` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
