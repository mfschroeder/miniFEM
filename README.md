# A minimal finite element program

In this repository we want to showcase the basic algorithm of discretizing a 
partial diferential equation (**PDE**) via the finite element method (**FEM**).
In our code we tried to stick as close as possible to the basic scheme of disrcetizing a PDE.
For the sake of simplicity we are only dealing with the most basic PDE, Possions problem, for this we
need to follow these steps.

![overview](https://user-images.githubusercontent.com/57176416/99834070-47b91280-2b63-11eb-8724-892ff9871833.png)

The attached code deals with the last step, the assembly of the resulting system of linear equations. We coverd all of these other steps on this [website](http://julianroth.org/documentation/fem/) where you can 
also find additional comments on the code.
For a general overview of FEM we would like to recommend our [YouTube Video](https://www.youtube.com/watch?v=P4lBRuY7pC4).
