## genepop

Le projet genepop a pour but la mise en place d'un package R qui permet l'utilisation
des fonctionnalités de l'outil [Genepop](http://kimura.univ-montp2.fr/~rousset/Genepop.htm) directement depuis l'environnement R.

Ainsi chacune des fonctionnalités fournies par l'outil Genepop est ici présente sous forme de fonction R.

L'outil Genepop étant développé en C++, la librairie Rcpp est utilisée afin de faire communiquer R et C++.

## Philosophie
L'idée derrière ce projet était d'identifier chacune des fonctionnalités de l'outil Genepop,
qui était utilisé en ligne de commandes ou à l'aide d'un menu interactif, et d'en extraire un ensemble
de fonctions ainsi que leurs signatures correspondantes.
Permettant ainsi à un utilisateur qui a installé notre package R, d'appeler directement ces fonctions depuis R.

## Architecture

```
├── genepop/
│   ├── man/
│   ├── inst/
│   │   ├── doc/
│   │   ├── genepop-shiny/
│   │   ├── extdata/
│   │   ├── extdoc/
│   │   ├── make-exe/
│   ├── R/
│   │   ├── Genepop.R
│   │   ├── RcppExports.R
│   ├── src/
│   │   ├── *.cpp
│   │   ├── *.h
│   ├── test/
│   │   ├── test-all.R
│   │   ├── testthat/
```

## Travail accompli
  - Réalisation des fonctions C++ qui représentent chaque fonctionnalité
  - Remplacement des appels à la fonction "exit" en C++ qui faisait planter R/RStudio
  - Initialisation/libération des variables globales
    - Afin d'initialiser et libérer correctement les variables globales pour utiliser plusieurs fonctionnalités successivement.
  - Mise à disposition d'un application web : genepop-shiny

## Dépendance
  Rcpp : Cette librairie permet d'utiliser les fonctions C++ depuis R

## Installation

    Installation depuis RStudio

   ```R
   library(devtools)
   library(git2r)

   devtools::install_git("http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop.git")

   #ou

   devtools::install_git("http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop.git", credentials = git2r::cred_user_pass ("your_login", "your_password"))
   ````


## Utilisation classique

Une fois que votre package genepop est installé.<br/>
Il suffit de choisir une fonctionnalité que l'on veut effectuer et de l'appeler avec les bons arguments.<br/>

Exemple pour un test de Hardy-Weinberg.<br/>
En entrée un [fichier](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/inst/extdata/sample.txt) au format Genepop<br/>
En sortie un fichier au format texte <br/>

  ```R
  #' @name Hardy-Weinberg
  #' @title Tests of Hardy-Weinberg genotypic proportions
  #' @description Compute variants of the exact conditional test for Hardy-Weinberg genotypic proportions.  The tests differ by their test statistics. \code{HWtable_analysis} handles a single table of genotype counts, and \code{test_HW} requires a standard genepop input file. See \href{../doc/all-menu-options.html#option-1-hardy-weinberg-hw-exact-tests}{this section} of the Genepop executable documentation for more information on the statistical methods.
  #' @param inputFile character: The path of the input file, in Genepop format
  #' @param which character: \code{'Proba'}, \code{'excess'}, and \code{'deficit'} to perform the probability test, score test for excess, and score tests for deficit, respectively, in each population and for each locus. \code{test_HW} additionally handles \code{'global excess'} and  \code{'global deficit'} for global tests for all loci and/or all populations, and \code{HWtable_analysis} additionally handles \code{'Fis'} to report basic information (allele frequencies and Fis).
  #' @param outputFile character: The path of the output file
  #' @param settingsFile character: The path of the settings file
  #' @param enumeration logical: whether to compute the complete enumeration test for samples with less than 5 alleles
  #' @param dememorization integer: length of dememorization step of Markov chain algorithm
  #' @param batches integer: Number of batches
  #' @param iterations integer: Iterations per batch
  #' @param verbose logical: whether to print some information
  #' @return The path of the output file is returned invisibly.
  #' @examples infile <- system.file('extdata', 'sample.txt',package='genepop')
  #' locinfile <- 'sample.txt' ## file in user's directory not in R's extdata directory
  #' check <- file.copy(infile,locinfile,overwrite=TRUE)
  #' test_HW(locinfile, which='deficit', 'sample.txt.D')
  test_HW(inputFile, which = "Proba", outputFile = "", enumeration = FALSE, dememorization = 10000, batches = 20, iterations = 5000, verbose = interactive()

  ```

## Utilisation avec fichier de paramètres

Certain fonction peuvent être utilisée avec un fichier de paramètres.<br/>

Une fois que votre package genepop est installé.<br/>
Il suffit de choisir une fonctionnalité que l'on veut effectuer avec un fichier de paramètres correcte.<br/>

Exemple pour un test de Hardy-Weinberg.<br/>
En entrée un [fichier](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/inst/extdata/sample.txt) au format Genepop et un fichier de [paramètres](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/inst/extdata/setting.txt)<br/>
En sortie un fichier au format texte <br/>

  ```R
  library(genepop)

  test_HW(inputFile, which = "Proba", outputFile = "", settingsFile = "", verbose = interactive())

  ```

## Sample

Le package genepop contient plusieurs [sample](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/tree/master/inst/extdata) qui vous aiderons à tester les différentes fonctionnalités si vous n'avez pas de fichier d'entrée conforme au format genepop.

  ```R
  #Permet de recupérer sous R, le chemin des différents fichier sample présent
  fpath <- system.file("extdata", "sample.txt", package="genepop")
  ```

## Documentation

Pour plus de détails concernent chaque fonction, vous référez à la [documentation](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/tree/master/man) R ou à la documentation de [genepop](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/tree/master/inst/doc)
<br/>

- **test_HW** : Tests of Hardy-Weinberg genotypic proportions ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/Hardy-Weinberg.Rd))

- **test_LD** : Tables and exact test for genotypic linkage disequilibrium ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/Linkage.Rd))

- **test_diff** : Tests of genic and genotypic differentiation ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/Differentiation.Rd))

- **struc** : Exact test on a single contingency table ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/Contingency-test.Rd))

- **Nm_private** : Private allele method ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/Nm_private.Rd))

- **basic_info** : Allele and genotype frequencies ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/basic_info.Rd))

- **genedivFis** : Gene diversities and Fis (or rho_IS) ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/genedivFis.Rd))

- **Fst** : Fst (or rho_ST) estimation ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/Fst.Rd))

- **ibd** : Isolation by distance ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/IBD.Rd))

- **conversion** : File conversions ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/conversion.Rd))

- **nulls** : Estimation of allele frequencies under genotyping failure. ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/nulls.Rd))

- **diploidize** : Various data manipulation utilities ([man](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/man/manipulation.Rd))


## Développement

- **Comment ajouter une fonctionnalité au package genepop** : <br/>
    Une fois que vous avez rajouté votre nouvelle fonctionnalité dans le code C++ de genepop, il va falloir la rendre disponible aux utilisateurs du package R. Pour cela il suffit de rajouter une fonction dans les fichiers [src/RGenepop.cpp](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/src/RGenepop.cpp) et [src/RGenepop.h](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/blob/master/src/RGenepop.h), en respectant le principe de fonctionnement des autres fonctions qui sont présentes , ainsi que la convention de nommage. <br/>

    L'étape suivante consiste à build votre package R. Une fois le build terminé, une fonction R est automatiquement créée dans le fichier R/RcppExports.R. Si votre fonction doit avoir des paramètres par défaut alors il faut la rajouter dans le fichier R/Genepop.R en lui rajoutant manuellement les valeurs par défauts de chaque paramètres. <br/> En respectant la aussi les conventions de nommage.

- **Utiliser genepop avec un fichier exécutable et non un package R** : <br/>
    Il suffit pour cela d'utiliser le fichier Makefile

## Shiny

Le package genepop propose aussi à ses utilisateurs un application web,  [genepop-shiny](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop/tree/master/inst/genepop-shiny), qui utilise [shiny](https://shiny.rstudio.com/) afin d'offrir une interface au package R
