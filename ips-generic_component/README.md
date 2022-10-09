
# generic_component.py
## 
*generic_component.py* is an IPS component wrapper for the case that the only thing needed
is the staging of input and state files to the work directory, execution of the physics
code, updating of state files and archiving of the output files.  It does no processing of
the input and state files nor processing of the output files or state files which have
been modified by the physics code, except possible copying of input or output files to
more convenient names. In other words it  is useful if the user has in hand, or other
components produce, input files readable by the physics code, and the code output files
can be used as is. It also supports running helper codes before and/or after running the
physics code. The component is generic in that it can, without modification, wrap any
physics code for use with IPS.  It could  also be used as a starting template for a more
general component.  The specification of  all aspects of the component to be implemented
is in the simulation configuration file.  For more details see the docstring in
*generic_component.py*.
