# JUNO

## Connect to remote computer

1.  Log to one of the multiple computer of the m2psa
	`ssh username@sbgli[1..20].in2p3.fr`
	- **N.B:** choose a number between 1 and 20. If an error occurs, try another number. In fact, many users can be connected, you have to choose a free one. You will have to enter your password later. Example:
	`ssh username@sbgli1.in2p3.fr`
1.  Log to a specific environment
	`ssh sbglmd[1..20]`
1.  Launch the ROOT environment
	`source /scratch/jcerasol/root_course/setup_root.sh`

## Clone a git repository
1.  Create a new repository wherever you want, let's call it `TIPP` 
	`mkdir TIPP`
1.  Init a git session
	`git init`
1.  Clone the git repository
	`git clone git@github.com:roronoarapha/JUNO.git`
	- **N.B:** Maybe you will have to generate a ssh key, so do :
	`ssh-keygen`

## Check you got this arborescence (at the minimum)

- JUNO
  - bin/
  - build/
  - include/
    - fonctions.h
	- plotter.h
  - src/
    - main.cpp
    - fonctions.cpp
	- plotter.cpp
  - output/
  - CMakeList.txt

## Check the branch
1. Go in **JUNO**, do `cd JUNO`
1. To see all the available branches, do :
`git branch -a`
1. The branch name with a star is the current branch. Example :

	> *main
	remotes/origin/HEAD -> origin/main
	remotes/origin/main
	remotes/origin/model_sigma
	remotes/origin/newfeature/rootgraph

1. (optional, your already are on the main branch) Got a remote branch in local. Example :
`git checkout -b newfeature/rootgraph origin/newfeature/rootgraph`

## Launch our code

1. Be sure you are in the **JUNO** folder
1. Go in the folder `build`
1. For a first use, trash its contain
	`rm -r *`
	- **N.B** Do not use this command line without thinking n times about what you want to do
1. From here, you can also trash the **output** folder, do
   `rm ../output/*`
1. Generate a Makefile associated to the project with
	`cmake ..`
1. Excute the code with 
	`make`

## This is the end

1. Check the output folder to got some graphics
1. Check your screen output to got some informations




