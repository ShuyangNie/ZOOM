.. _usage:

----------
User Guide
----------

``zoom`` dedicates to provide a Python platform for integrating multimodal, multiscale neuro‑omics data types. 
In our present work, we have successfully linked spatial brain phenotypes (SBPs), the anatomically comprehensive 
Allen Human Brain Atlas (AHBA), and single‑cell RNA sequencing (scRNA-seq) datasets.

This user guide describes, in sequence, how to preprocess the above data types with ``zoom``, implement the 
imaging‑transcriptomics paradigm which associate SBPs with AHBA, and finally compute single‑cell enrichment 
scores for SBP‑relevant gene sets. If you have further questions, please consult :ref:`api_ref` or submit 
an issue on `GitHub <https://github.com/SpaTrek/ZOOM/issues>`_.

.. toctree::
   :caption: Table of Contents
   :maxdepth: 2

   user_guide/prepare.rst
   user_guide/imt.rst
   user_guide/sc_score.rst
