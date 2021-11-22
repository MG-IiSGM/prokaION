# prokaion

Pipeline for analysis of any Prokaryote MinION samples

For the development of this pipeline it is necessary to download [Guppy](https://community.nanoporetech.com/downloads/guppy/release_notes) from [Nanopore](https://community.nanoporetech.com/downloads) for the basecalling and de-multiplexing steps. If you have used [MinKNOW](https://community.nanoporetech.com/downloads/minion_release/release_notes) to obtain the .fast5 files, you can also perform these steps within the same software.


To install the environment with all the specifications installed, proceed as follows:
  conda create -n 'env_name'
  conda install --yes --file 'spec_file' -c conda-forge
