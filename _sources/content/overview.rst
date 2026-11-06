.. _overview:

Overview
========

The CATchment HYdrology (CATHY) model couples a finite element solver for the Richards equation describing ow in variably saturated porous media and a finite difference solver for the diffusion wave equation describing surface flow propagation throughout a hill slope and stream channel network identified using terrain topography and the hydraulic geometry concept 
The mathematical model is described by a system of two differential equations :cite:p:`Camporese2010`:

1. Subsurface flow equation:

.. math:: 
	
   S_{w}S_{s}\frac{\partial \psi}{\partial t} \,+\,\Phi\frac{\partial s_{w}}{\partial t} =\,\nabla.\,\left[K_{s}K_{r}\left(\nabla\psi\,+\,n_{z}\right)\right]\,+\,q_{ss}


2. Surface flow equation:

.. math::
	
   \frac{\partial Q}{\partial t} + c_{k}\frac{\partial Q}{\partial s} = D_{h}\frac{\partial^{2} Q}{\partial s^{2}} + c_{k}q_{s}


The 3-D Richards equation is solved numerically by Galerkin finite elements in space using tetrahedral elements and linear basis functions, and by a weighted finite difference scheme for integration in time. The nonlinear characteristic relationships :math:`K_{r}(\Psi)` and  :math:`S_{w}(\Psi)` are specified using either van Genuchten and Nielsen [1985], Brooks and Corey [1964], or Huyakorn et al. [1984] expression. Linearization via Newton or Picard iteration is used in the solution procedure. 

