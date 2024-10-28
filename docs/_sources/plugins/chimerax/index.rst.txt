*************
Find Cavities
*************

The **pyKVFinder** package is integrated into `UCSF ChimeraX <https://www.cgl.ucsf.edu/chimerax/>`_ as the **Find Cavities** tool, available from version 1.9. This tool allows the user to run pyKVFinder in ChimeraX and visualize the results directly on the ChimeraX viewer.

Tool
====

The **Find Cavities** tool can be accessed from:

* ``Tools`` > ``Binding Analysis`` > ``Find Cavities`` 
  
or,

* ``Tools`` > ``Structure Analysis`` > ``Find Cavities``

For further details on Find Cavities tool, please refer to https://www.cgl.ucsf.edu/chimerax/docs/user/tools/findcavities.html.

Command
=======

The Find Cavities tool uses the **kvfinder** command. The command can be used in the ChimeraX command line interface (CLI) as follows:

.. code-block::

    kvfinder  model-spec  [ probeIn  r1 ] [ probeOut  r2 ] [ removalDistance  d ] [ volumeCutoff  minvol ] [ gridSpacing  s ] [ showTool  true | false ] [ surfaceType  SAS | SES ] [ boxOrigin  x,y,z ] [ boxExtent  length | lx,ly,lz ]

For further details on kvfinder command, please refer to https://www.cgl.ucsf.edu/chimerax/docs/user/commands/kvfinder.html.
