.. _api_ref:

.. currentmodule:: zoom

API Documents
=============

.. contents:: **List of modules**
   :local:

.. _ref_core:

:mod:`zoom.core` - Core Class
-----------------------------
.. automodule:: zoom.core
   :no-members:
   :no-inherited-members:

.. currentmodule:: zoom

.. autosummary::
   :template: class.rst
   :toctree: generated/

   zoom.core.ZOOM
   zoom.core.ZOOM_SC

.. _ref_data:

:mod:`zoom.data_loader` - Data Loader
-------------------------------------
.. automodule:: zoom.data_loader
   :no-members:
   :no-inherited-members:

.. currentmodule:: zoom

.. autosummary::
   :template: function.rst
   :toctree: generated/

   zoom.data_loader.load_sc
   zoom.data_loader.load_gii
   zoom.data_loader.load_parc

.. _ref_prepare:

:mod:`zoom.prepare` - Prepare Input Data
----------------------------------------
.. automodule:: zoom.prepare
   :no-members:
   :no-inherited-members:

.. currentmodule:: zoom

.. autosummary::
   :template: function.rst
   :toctree: generated/

   zoom.prepare.abagen_ctx
   zoom.prepare.process_SBP
   zoom.prepare.fetch_medial_wall

.. _ref_pls:

:mod:`zoom.pls_tool` - Imaging Transcriptomics
----------------------------------------------
.. automodule:: zoom.pls_tool
   :no-members:
   :no-inherited-members:

.. currentmodule:: zoom

.. autosummary::
   :template: function.rst
   :toctree: generated/

   zoom.pls_tool.optimal_component_eval
   zoom.pls_tool.run_component_eval
   zoom.pls_tool.model_eval
   zoom.pls_tool.pls_perm
   zoom.pls_tool.vip_perm
   zoom.pls_tool.boot_pls1
   zoom.pls_tool.pls1_perm

.. _ref_sc:

:mod:`zoom.sc_tool` - Single Cell Scoring
-----------------------------------------
.. automodule:: zoom.sc_tool
   :no-members:
   :no-inherited-members:

.. currentmodule:: zoom

.. autosummary::
   :template: function.rst
   :toctree: generated/

   zoom.sc_tool.preprocess
   zoom.sc_tool.rank_expression
   zoom.sc_tool.compute_gss
   zoom.sc_tool.select_ctrl
   zoom.sc_tool.score_cell_zoom
   zoom.sc_tool.group_bh
   zoom.sc_tool.downstream_DEG
   zoom.sc_tool.downstream_region_enrich
   zoom.sc_tool.gsea_perm