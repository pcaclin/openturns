#! /usr/bin/env python

from __future__ import print_function
import openturns as ot


eps = 1e-2
# Instance creation
myFunc = ot.Function(['x1', 'x2'], ['f1', 'f2', 'f3'], [
                                  'x1*sin(x2)', 'cos(x1+x2)', '(x2+1)*exp(x1-2*x2)'])
print("myFunc (before substitution) = ", repr(myFunc))
epsilon = ot.NumericalPoint(myFunc.getInputDimension(), eps)
inPoint = ot.NumericalPoint(epsilon.getDimension(), 1.0)
myGradient = ot.CenteredFiniteDifferenceGradient(
    epsilon, myFunc.getEvaluation())

print("myGradient=", repr(myGradient))
print("myFunc.gradient(", repr(inPoint), ")=",
      repr(myFunc.gradient(inPoint)))
print("myGradient.gradient(", repr(inPoint), ")=",
      repr(myGradient.gradient(inPoint)))

# Substitute the gradient
myFunc.setGradient(myGradient)
print("myFunc (after substitution) = ", repr(myFunc))

print("myFunc.gradient(", repr(inPoint), ")=", repr(
    myFunc.gradient(inPoint)), " (after substitution)")
