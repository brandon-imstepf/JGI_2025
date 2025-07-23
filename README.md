# BBTools-ML: Machine Learning Extensions for Gene Prediction

This is a fork of Brian Bushnell's [BBTools suite](https://jgi.doe.gov/data-and-tools/bbtools/), developed over eight weeks in Summer 2025. This project adds a machine learning pipeline for prokaryotic gene prediction.

The primary additions are:
1.  Tools for processing `.gff` files into a machine-learning-readable format (`.tsv`).
2.  A custom **Convolutional Neural Network (CNN)**, written entirely in Java with no external ML packages.
3.  Integration of trained models into `callgenes.sh` to filter and improve gene predictions.

---
## Key Features

* **Data Generation**: Convert annotated genomes (FASTA + GFF) into feature vectors for training.
* **Custom CNN Trainer**: Train a convolutional neural network directly within the Java environment.
* **Inference Mode**: Use a pre-trained `.bbnet` model within `callgenes.sh` to score and filter gene candidates, improving prediction accuracy.

---
## Setup

This project is built using the standard BBTools files From the main `bbmap/` directory, run:

```bash
./compile.sh
```
If you ever add more files to any of the directories, add them to `sources.txt` as `compile.sh` grabs the files to build from there.

## Example Workflow

The process involves three main steps: generating a training set, training a model, and running inference.

### Step 1: Generate a Training Set

Use the `gff2tsv.sh` script to convert a GFF file and its corresponding FASTA file into a `.tsv` feature vector file. The `buildtrainingset.sh` script automates this for multiple files.


# Example for a single genome
```bash
gff2tsv.sh in=genome.fna.gz gff=known_genes.gff.gz out=training_data.tsv
```

# Example using the helper script for a folder of genomes
```bash
buildtrainingsetfolder.sh in=genomes/ out=training_set.tsv
```

### Step 2: Train a Network

You can use the generated `.tsv` file to train a model. This fork provides a built-in CNN trainer.

# Train the custom Java CNN
```bash
cnntrain.sh in=training_set.tsv out=my_model.bbnet epochs=10
```
### Step 3: Run Inference with CallGenes

Use the modified `callgenes.sh` with a trained model (`.bbnet` file) to predict genes on a new genome.

```bash
callgenes.sh in=new_genome.fna out=predictions.gff net=my_model.bbnet cutoff=0.8
```
* net=<file>: Activates the neural network filter.
* cutoff=<float>: Sets the score threshold for keeping a gene prediction.

## File Descriptions

### New Java Classes

* `ml/CNNNetwork.java`, `ml/CNNTrainer.java`: Core classes for the custom CNN.
* `ml/Layer.java` & subclasses: Defines the layers used in the network (`ConvolutionLayer`, `MaxPoolingLayer`, `DenseLayer`, `OutputLayer`).
* `ml/BinaryCrossEntropyLoss.java`: The loss function for training.
* `prok/CallGenesHelper.java`: Helper methods for `CallGenes`, including feature generation.
* `prok/GfftoTSV.java`: The underlying program for converting GFF to TSV format.

### New Shell Scripts

* `buildtrainingset.sh`, `buildtraininsetfolder.sh`: Scripts to automate the creation of training data.
* `gff2tsv.sh`: Wrapper for `GfftoTSV`.
* `gffsetop.sh`: Helper script for GFF file operations.
* `cnntrain.sh`: Wrapper script for the `CNNTrainer`.
* `callgenesfolder.sh`: Helper script to run `callgenes.sh` on a folder of genomes.

### Modified Files

* `prok/CallGenes.java`: Modified to support training data generation and inference with a neural network.
* `callgenes.sh`: Usage message updated with new parameters.

