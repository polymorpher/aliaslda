### Created and maintained by: [Aaron Li](https://www.linkedin.com/in/aaronqli)

A demonstration of AliasLDA, AliasPDP, and AliasHDP, as described in the KDD 2014 paper [Reducing the Sampling Complexity of Topic Models](https://storage.googleapis.com/aliaslda/kdd2014afterPeerReview.pdf) (Aaron Li, Amr Ahmed, Sujith Ravi, Alexander Smola). 

Slides (with presentation captions), presented at KDD 2014 conference, is available [here](https://storage.googleapis.com/aliaslda/kdd2014talk.pdf). The slide themselves are kindly created by Alex Smola two days before the presentation, to help me focus on the presentation process itself.

The source code is identical to what was used in the original experiments, except for 

1. a CMake + CLion project wrapper
2. a small patch on SparseLDA, to fix a bug.

With the bug, topics with documented-count 1 are sampled slightly incorrectly. 

In the paper, the effect is visible in the experiment results: SparseLDA has slightly higher perplexity after convergence than AliasLDA. 

However, it can be verified that the sampling speed of SparseLDA are the same before and after this patch, therefore all results in the paper still hold. 

Q&A
=========
Q: The code looks really messy and hacky! Can you clean up the code?

A: That's because I wrote it very quickly (2-3 weeks in total) for doing experiments, and never worked on it again since Feb 2012. I will try to clean it up when I have time.

-----

Q: Are you planning to make a library package available?

A: I would love to! I am currently spending 100+ hours per week working on [**Scaled Inference**](https://scaledinference.com/). But I will make time when I can

-----

Q: How about the latest awesome techniques discovered and published by other people? e.g: [LightLDA](http://arxiv.org/abs/1412.1576), [WarpLDA](http://arxiv.org/abs/1510.08628), [F-tree](https://en.wikipedia.org/wiki/Fenwick_tree?oldformat=true)-based-LDA, etc.?

A: They are great. I learned a lot from these papers and I can see techniques that I can use to make this package much better. I will integrate these techniques and update this package later.

-----

Q: What about large-scale version, such as [High Performance Latent Variable Models](http://arxiv.org/abs/1510.06143v1)?

A: Although I designed the large scale version, wrote all the code and did all experiments behind that paper, I do not own the code. Google does. So unfortunately it won't be available here.

-----

Q: How about you writing a version using one of the open source large-scale machine learning framework, such as [Parameter Server](https://github.com/dmlc/ps-lite), or [Petuum](https://petuum.github.io/index.html)?

A: Great idea! I will do that when I have time.

Instructions
=========
I recommend using an IDE to compile this, such as [CLion](https://www.jetbrains.com/clion/). You can use the run-configuration I provided in the project, to start quickly.

The default CMakeLists.txt only compiles AliasLDA and SparseLDA related stuff. You can uncomment some lines to compile other implementations, such as HDPLDA, PDPLDA, AliasPDP, AliasHDP, and some hacky tests, etc.

The **main** function is in **AliasLDATester.cpp**, you can see the command-line parameters there, or play with it in whatever way you want.

Example (assuming the binary is compiled as `aliaslda`):

```!sh
aliaslda ~/aliasdata/conf/test.ext1.conf aliaslda 0.1 0.1 1024 100
```

This launches an experiment using the config file `~/aliasdata/conf/test.ext1.conf`, with alpha=0.1, beta=0.1, 1024topics, and 100 iterations.

`~/aliasdata/conf/test.ext1.conf` is the config file describing the data locations. See next section for detail.

Data
========
I shared the data from [my Dropbox folder](https://www.dropbox.com/sh/zryf092lcatwtc0/AABybC3JS7pAg27LOcd6TawKa?dl=0). In `test` folder, you can find the config file mentioned in last section. In `data` folder, you can find actual data.

To use them directly without any modification, just copy them to the root folder of this project.

The config file is pretty self-explanatory. The paths defined in the config are relative to run-time path (i.e your $PWD when you run the binary)

There are two formats for the data files (as you can see in config files).

- Standard: Every line is a document. Every token (word) in the line is an integer (in text form). The mapping between the integer token and the text is provided in the vocabulary file. 
- Bag: The three lines at the beginning are meaningless. After that, each line is in the form of `<docId> <tokenId> <counts>`

 




