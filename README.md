# A minimal finite element program

In this repository we want to showcase the basic algorithm of discretizing a 
partial diferential equation (**PDE**) via the finite element method (**FEM**).
In our code we tried to stick as close as possible to the basic scheme of discretizing a PDE.
For the sake of simplicity, we are only dealing with the most basic PDE, Possion's problem. For this we
need to follow these steps:

![overview](https://user-images.githubusercontent.com/57176416/99834070-47b91280-2b63-11eb-8724-892ff9871833.png)

The code in this repository deals with the last step, the assembly of the resulting system of linear equations. We explained all of the other steps on this [website](http://julianroth.org/documentation/fem/), where you can also find additional comments on the code.
For a general overview of FEM we would like to recommend our [YouTube Video](https://www.youtube.com/watch?v=P4lBRuY7pC4).

## Numerical experiments

### Problem statement
Find <img src="https://latex.codecogs.com/gif.latex?u:&space;\Omega&space;\rightarrow&space;\mathbb{R}" title="u: \Omega \rightarrow \mathbb{R}" /> such that <br><br>
<img src="https://latex.codecogs.com/gif.latex?-\Delta&space;u&space;=&space;-1&space;\text{&space;in&space;}&space;\Omega&space;:=&space;(0,1)&space;\times&space;(0,1)&space;\\&space;\phantom{aaaaaa}&space;u&space;=&space;0&space;\text{&space;on&space;}&space;\partial\Omega" title="-\Delta u = -1 \text{ in } \Omega := (0,1) \times (0,1) \\ \phantom{aaaaaa} u = 0 \text{ on } \partial\Omega" />

### Numerical solution
Discretizing with linear finite elements and mesh size h = 0.05, we get this FEM solution:

![solution](https://user-images.githubusercontent.com/42407091/99839508-19d7cc00-2b6b-11eb-8b39-ba331b80b285.png)
