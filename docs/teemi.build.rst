Build module
============

The primary actively maintained build APIs are the PCR and transformation helpers.
The robot assembly and picklist modules below are kept for backwards compatibility
with older Flowbot-oriented workflows and legacy notebooks.

PCR module
----------

.. automodule:: teemi.build.PCR
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:



Legacy robot assembly module
----------------------------

.. note::

   ``teemi.build.robot_assembly`` is a legacy compatibility surface.
   New code should prefer newer build-planning workflows rather than relying on
   the old robot/picklist abstractions.

.. automodule:: teemi.build.robot_assembly
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Transformation module
---------------------

.. automodule:: teemi.build.transformation
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Legacy containers and picklists module
--------------------------------------

.. note::

   ``teemi.build.containers_wells_picklists`` is a legacy compatibility surface.
   It remains available for historical workflows, tests, and older notebooks.

.. automodule:: teemi.build.containers_wells_picklists
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
