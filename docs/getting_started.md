# Getting Started

This section contains a quick guide for installing, and using basespaceapp. For the  detailed syntax guide, see the [Basespace Syntax][syntax]

[syntax]: syntax.html


## Installation

Basespaceapp is available on PyPi so installing is easy.

```bash 
$ pip install basespaceapp
```

 
To make it easy, basespaceapp runs on Python 2.7, 3.4, and 3.5!


## Using basespaceapp

By default, basespaceapp is both a command line tool and a Python module that can be used to build and run pipelines in code. This means that whether you're a user, or a developer Basespaceapp can be adapted to fit your needs.

To see basespaceapp's help menu, type the following, just as you'd expect.

```bash 
$ basespaceapp --help
```

## Sample Pipeline

Here's a simple pipeline you can use for testing basespaceapp. Typically complex pipelines are used for things like bioinformatics, or some batch processing

But first, we need some sample files to work with. Run these commands to generate them.

```bash 
$ echo "SAMPLE DATA 1" > test_file.1.txt
$ echo "SAMPLE DATA 2" > test_file.2.txt
$ echo "SAMPLE DATA 3" > test_file.3.txt
```

Now that we have our data, let's analyze it! Here's our sample pipeline:

```bash
[COMMANDS]
# Remove the ending number from each of our data files.
cut -f 1-2 -d ' ' {1||2||3} > {o}

# Paste each of the files together and save it to a final output.
# Since this is our last step, and only 1 output there's no need to have 
# basespaceapp name the output file. We'll call it something ourselves. 
paste {1.1,1.2,1.3} > final_output.txt

[FILES]
1. test_file.1.txt
2. test_file.2.txt
3. test_file.3.txt
```

Save that as `sample_pipeline.mp`, open a terminal, and `cd` to that directory.


### Run the sample pipeline locally

Local execution is the default for basespaceapp so you just need to specify your basespaceapp file and an output destination.

```bash
$ basespaceapp -o pipeline.sh sample_pipeline.mp
```

This will generate an output script named `pipeline.sh` which will run the pipeline. Simply run it to start your pipeline!

```bash
$ sh pipeline.sh
```

That's it! Basespaceapp will run in the foreground watching your jobs complete until everything finishes.


### Run the sample pipeline on PBS

Simply change the basespaceapp command to the following:

```bash
$ basespaceapp -o pipeline.sh -j pbs sample_pipeline.mp
```

Then simply submit basespaceapp as a job:

```bash
$ qsub pipeline.sh
```

Basespaceapp will run as a job on the PBS/Torque queue and submit other jobs to the same queue! It will keep tabs on the running jobs and submit them when they're ready, then exit when all jobs finish.
