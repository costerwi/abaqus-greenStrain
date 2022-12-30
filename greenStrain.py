"""Abaqus script to compute Green strain from Nominal strain

Green's strain is sometimes used in large-rotation, small-strain formulations
for such problems as shell buckling. It is not available directly as an Abaqus
field output nor is the displacement gradient which is commonly used in its
calculation.

This script computes and stores Green strain (GE) field output for all odb step
frames which include Nominal strain (NE).

Reference:

    https://www.continuummechanics.org/greenstrain.html
    https://help.3ds.com/2023/english/dssimulia_established/SIMACAETHERefMap/simathe-c-strainmeas.htm

Execution:

    abaqus python greenStrain.py Job-1.odb [Job-2.odb ...]

Optional automatic execution:

    Copy the onJobCompletion() method into a local abaqus_v6.env file for
    automatic execution whenever an Abaqus job completes. Make sure this
    script is in the working directory or the PYTHONPATH

Green strain is calculated based on the relationships:

    GE = 1/2(F*FT - I)
    F*FT - I = V*V - I
    V = NE + I

    Where
    GE is the Green strain tensor in current configuration
    F is the displacement gradient tensor
    V is the stretch tensor in current configuration
    NE is the Nominal strain tensor in current configuration

Carl Osterwisch, December 2022
https://github.com/costerwi/abaqus-greenStrain
"""

from __future__ import print_function, with_statement
import numpy as np

__version__ = "1.0.0"
enableLocalCoordSystem = True  # set False to improve performance if not needed


def onJobCompletion():
    """Copy this method into abaqus_v6.env for automatic execution"""
    import os, greenStrain
    greenStrain.fromOdb(os.path.join(savedir, id + ".odb"))


def quat2matrix(quatArray):
    """Convert array of quaternions into array of rotation matrices
 
    >>> A = quat2matrix([[0.0, 0.0, 0.0, 1.0],
    ...              [0, 0, np.sin(np.pi/4), np.cos(np.pi/4)]])
    >>> print(A.shape)
    (2, 3, 3)
    >>> A
    array([[[ 1.,  0.,  0.],
            [ 0.,  1.,  0.],
            [ 0.,  0.,  1.]],
    <BLANKLINE>
           [[ 0., -1.,  0.],
            [ 1.,  0.,  0.],
            [ 0.,  0.,  1.]]])
    """

    quatArray = np.asarray(quatArray)
    assert 2 == len(quatArray.shape), "Must be a 2D array"
    assert 4 == quatArray.shape[1], "Must be an array of quaternions"
    q1, q2, q3, q0 = quatArray.T
    return np.transpose([
        1 - 2*q2*q2 - 2*q3*q3,  2*q1*q2 - 2*q0*q3,  2*q1*q3 + 2*q0*q2,
        2*q1*q2 + 2*q0*q3,  1 - 2*q1*q1 - 2*q3*q3,  2*q2*q3 - 2*q0*q1,
        2*q1*q3 - 2*q0*q2,  2*q2*q3 + 2*q0*q1,  1 - 2*q1*q1 - 2*q2*q2
    ]).reshape(-1, 3, 3)


def tensor2(A):
    """Compute dot product of array of tensor results A with itself

    >>> E11, E22, E33, E12, E13, E23 = 0.1, 0.2, 0.3, 0.4, 0.5, 0.6
    >>> tensor2( [[E11, E22, E33, E12, E13, E23],
    ...           [0.2, 0.0, 0.0, 0.3, -.5, 0.0]] )
    array([[ 0.42,  0.56,  0.7 ,  0.42,  0.44,  0.5 ],
           [ 0.38,  0.09,  0.25,  0.06, -0.1 , -0.15]])
    """

    A = np.asarray(A)
    assert 2 == len(A.shape), "Data must be a 2D array of tensors"
    assert 6 == A.shape[1], "Data elements must be in Abaqus symmetric tensor format"
    result = np.empty_like(A)
    t11, t22, t33, t12, t13, t23 = range(6)  # column indices for tensor array
    t21, t31, t32 = t12, t13, t23  # symmetry

    result[:,t11] = A[:,t11]*A[:,t11] + A[:,t12]*A[:,t21] + A[:,t13]*A[:,t31]
    result[:,t22] = A[:,t21]*A[:,t12] + A[:,t22]*A[:,t22] + A[:,t23]*A[:,t32]
    result[:,t33] = A[:,t31]*A[:,t13] + A[:,t32]*A[:,t23] + A[:,t33]*A[:,t33]
    result[:,t12] = A[:,t11]*A[:,t12] + A[:,t12]*A[:,t22] + A[:,t13]*A[:,t32]
    result[:,t13] = A[:,t11]*A[:,t13] + A[:,t12]*A[:,t23] + A[:,t13]*A[:,t33]
    result[:,t23] = A[:,t21]*A[:,t13] + A[:,t22]*A[:,t23] + A[:,t23]*A[:,t33]
    return result


def green(NEarray):
    """Calculate Green strain array of tensors from NE array of tensors

    >>> NE11, NE22, NE33, NE12, NE13, NE23 = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0
    >>> green( [[NE11, NE22, NE33, NE12, NE13, NE23],
    ...         [ 0.2,  0.0,  0.0,  0.3,  -.5,  0.0]] )
    array([[ 6.625, 10.5  , 15.125, 17.5  , 21.   , 26.   ],
           [ 0.262,  0.011,  0.031,  0.33 , -0.55 , -0.037]])
    """

    NEarray = np.asarray(NEarray)
    assert 2 == len(NEarray.shape), "NEarray must be a 2D array of tensors"
    assert 6 == NEarray.shape[1], "NEarray elements must be in Abaqus symmetric tensor format"
    V = NEarray.copy()
    # V = NE + I
    V[:, 3:] /= 2  # NE engineering shear to tensor shear
    V[:, :3] += 1  # add identity matrix
    # GE = 1/2(V*V - I)
    V = tensor2(V)  # compute V*V
    V[:, :3] -= 1  # subtract identity matrix
    V[:, :3] /= 2  # Green strain with shear in engineering form
    return V


def calculateGE(NE, outputFrame):
    """Calculate Green strain (GE) fieldOutput and store in outputFrame"""

    GE = outputFrame.FieldOutput(
        name="GE",
        description="Green strain components",
        type=NE.type,
        validInvariants=NE.validInvariants,
        isEngineeringTensor=True,
    )

    for block in NE.bulkDataBlocks:
        options = dict(
            position=block.position,
            instance=block.instance,
            labels=block.elementLabels,
            data=green(block.data),
        )
        if np.any(block.sectionPoint):
            options["sectionPoint"] = block.sectionPoint
        if enableLocalCoordSystem and np.any(block.localCoordSystem):
            # Abaqus 2023 returns but does not accept quaternion localCoordSystem.
            # The quaternion array must be converted to a matrix array.
            quatArray = block.localCoordSystem.copy()
            quatArray[:, 3] *= -1  # Inverse
            options["localCoordSystem"] = np.ascontiguousarray(quat2matrix(quatArray))
        GE.addData(**options)


def fromOdb(odbName):
    """Add Green strain to each odb frame which contains Nominal strain"""

    from odbAccess import openOdb
    from contextlib import closing

    with closing(openOdb(odbName)) as odb:
        for step in odb.steps.values():
            for frame in step.frames:
                if "NE" in frame.fieldOutputs:
                    NE = frame.fieldOutputs["NE"]
                    calculateGE(NE, frame)


import sys
for arg in sys.argv[1:]:
    if "--help" == arg:
        print(__doc__)
    elif "--test" == arg:
        import doctest
        np.set_printoptions(precision=3, suppress=True)
        doctest.testmod(verbose=True)
    else:
        print(arg, file=sys.stderr)
        fromOdb(odbName = arg)
