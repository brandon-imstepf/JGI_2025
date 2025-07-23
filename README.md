A fork of Brian Bushnell's BBTools, made with effort over eight weeks in Summer 2025 to create tools for processing .gff files into a Machine-Learning-readable format (.tsv), then train a Feed Forward Neural Network on it. Finally, I developed a Convolutional Neural Network entirely within Java with no external packages to also train on the same data.

Added Files:
ml.CNNNetwork
ml.CNNTrainer
ml.Layer
ml.ConvolutionLayer
ml.MaxPoolingLayer
ml.DenseLayer
ml.OutputLayer
ml.BinaryCrossEntropyLoss

prok.CallGenes*
prok.CallGenesHelper
prok.GfftoTSV

buildtrainingset.sh
buildtraininsetfolder.sh
gff2tsv.sh
gffsetop.sh
cnntrain.sh
callgenes.sh*
callgenesfolder.sh


(* denote edited files, not created by me.)
