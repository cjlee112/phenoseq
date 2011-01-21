-------------
Steps
-------------

* Demixing
* Aligning
* Analyzing

---------------
Demixing
---------------

The first task is to demultiplex the raw pooled sequence fragments. 
Let's assume the data files have the format ``s_L_{1,2}_n_qseq.txt``. 
The number L is the lane, and the number n is the tile number. 
For each tile, there are 2 files. The fragments are in the file with a 
1 after the second underscore and the tags associated to those fragments 
are in the file with the 2 after the underscore. The tags used can be 
found in the file "SampleSheet.csv". Because the tags are designed with 
error correction in mind, if a reported tag differs from one of the tags 
in the sample sheet file by just one base, we assume that it was tagged 
with that tag.

The script demultiplex.py can process the raw files and separate the reads 
by tag::

    python demultiplex.py SampleSheet.csv data_directory

The output is a collection of files of the form ``s_(lane)_(tag)_qseq.txt`` .

----------------
Aligning
----------------

The resulting demixed reads can be aligned with novoalign. First we must build an indexed version of a refence genome. If the reference genome file is named ``NC_000913.fna``, the command to build the index is::

    ./novocraft/novoindex ecoli.nix NC_000913.fna

which creates the index ``ecoli.nix``. 

Alignment then proceeds with the command::

    ./novoalign/novocraft/novoalign -d ecoli.nix -f s_1_ACTTGA_qseq.txt

which aligns the reads in the file ``s_1_ACTTGA_qseq.txt`` against the index. Since there will be several such files, we can automate this process with a simple shell script (align_batch.csh)::

    #!/bin/csh 
    #                                                                                                          
    foreach l ( 1 2 3 )                                                                                        
    foreach t ( ACAGTG ACTTGA ATCACG CAGATC CGATGT CTTGTA GATCAG GCCAAT TGACCA TTAGGC )                      
        ./novoalign/novocraft/novoalign -d ecoli.nix -o SAM -f s_$l"_"$t"_"qseq.txt > aligned_s_$l"_"$t.sam    
    end                                                                                                      
    end

The ``-o SAM`` option outputs the data in SAMTOOLS format.

See 

http://www.novocraft.com/wiki/tiki-index.php?page=Getting+Started&structure=Novocraft+Technologies&page_ref_id=70

for more information on using novoalign.

-----------------
Analyzing
-----------------

SNP analysis is handled by mkpileup.py and analyze.py. The first script produces vcf files from the .sam files from the alignment::

    python mkpileup.py NC_000913.gbk [ACGT]*.vcf. 

These files can then be processed by analyze.py::

    python analyze.py NC_000913.gbk [ACGT]*.vcf

