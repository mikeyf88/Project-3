#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:34:19 2021

@author: kendrick
"""

import numpy as np

# compute unknown displacements 
def ComputeDisplacements(K, F, n_unknowns):
    # extract submatrix of unknowns
    K11 = K[0:n_unknowns,0:n_unknowns]
    F1 = F[0:n_unknowns]
    
    d = np.linalg.solve(K11,F1)
    
    return d

# postprocess the forces at known displacement nodes
def PostprocessReactions(K, d, F, n_unknowns, nodes):
    # These are computed net forces and do not
    # take into account external loads applied
    # at these nodes
    F = np.matmul(K[n_unknowns:,0:n_unknowns], d)
    
    # Postprocess the reactions
    for node in nodes:
        if node.xidx >= n_unknowns:
            node.AddReactionXForce(F[node.xidx-n_unknowns][0] - node.xforce_external)
        if node.yidx >= n_unknowns:
            node.AddReactionYForce(F[node.yidx-n_unknowns][0] - node.yforce_external)
        
    return F

# determine internal member load
def ComputeMemberForces(bars):
    for bar in bars:

        E = bar.E
        A = bar.A
        L = bar.Length() * 12  # Convert feet to inches

        [lambdax, lambday] = bar.LambdaTerms()

        init_node = bar.init_node
        end_node = bar.end_node

        u1 = init_node.xdisp
        v1 = init_node.ydisp
        u2 = end_node.xdisp
        v2 = end_node.ydisp

        axial_load = (E * A / L) * (
            -lambdax * u1
            -lambday * v1
            +lambdax * u2
            +lambday * v2
        )

        bar.axial_load = axial_load

    return bars
   
        
# compute the normal stresses
def ComputeNormalStresses(bars):
    # Compute normal stress for all bars

    for bar in bars:

        # Extract axial force and cross-sectional area
        axial_load = bar.axial_load
        A = bar.A

        # Calculate normal stress
        normal_stress = axial_load / A

        # Store normal stress in bar
        bar.normal_stress = normal_stress

    return bars

# compute the critical buckling load of a member
def ComputeBucklingLoad(bars):
    for bar in bars:

        E = bar.E
        L = bar.Length() * 12  # Convert feet to inches

        K = 1.0

        I_strong = bar.GetStrongSecondMoment()
        I_weak = bar.GetWeakSecondMoment()

        I = min(I_strong, I_weak)

        Pcr = (np.pi**2 * E * I) / ((K * L)**2)

        bar.buckling_load = Pcr

    return bars

    return bars
