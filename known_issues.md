# Known issues pyCATHY DA

- ERT assimilate in non parallel model lead to wrong observation size (missing a stack somewhere)


# Known issues pyCATHY



- FP arg to update soil file is a dict or a pd? as for now a dict but we should be able to pass directly a pandas dataframe

  - just pass df.to_dict(orient='list')

  

## Parm file

- non-linear convergence induce change in delta T (visible on risul file)


## while running

the processor only return :  nsf  (# of seepage faces)               =      0

--> problem with the mesh
- Discretisation of the mesh in the z direction +1



Program received signal SIGSEGV: Segmentation fault - invalid memory reference.
Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG
--> bad atmbc float number



root depth is not an element of the mesh layer!


 ELEMENTO DIAGONALE DI L NULLO,I,J =    1      -0.32019E-08

 
