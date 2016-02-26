# Basespaceapp

*A template to make Basespace native apps.*

Basespaceapp is a simple template for building and running Basespace native apps. If you use a PBS/Torque queue for cluster computing, or if you have complex batch processing that you want simplified, metapipe is the tool for you.

Basespaceapp's goal is to improve **readability**, and **maintainability** when building complex Basespace native apps.

In addition to helping you generate and maintain complex Basespace native apps, **basespaceapp also helps you debug them**! How? Well basespaceapp watches your app execute and keeps tabs on them. This means, unlike conventional batch queue systems like PBS/Torque alone, metapipe can give you accurate error information, and even resubmit failing jobs! Basespaceapp enhances the power of any PBS/Torque queue! 

- What if I [don't use PBS/Torque](#other-queue-systems), or [a queue system at all?](#no-queue-no-problem)


## What does it do?

In the bad old days (before basespaceapp), if you wanted to make a Basespace native apps, you needed to know how to code json input form specifications and javascript callback functions. **Not anymore!** Basespaceapp makes it easy to build and run your analysis pipelines! **No more code, just commands!** This makes your pipelines easy to understand and change! 

A sample basespace file can be found in [Basespaceapp Syntax](syntax.html)


## No Queue? No Problem!

Lots of people don't use a PBS/Torque queue system, or a queue system at all, and basespaceapp can help them as well! Basespaceapp runs locally and will give you all the same benefits of a batch queue system! It runs jobs in parallel, and provide detailed feedback when jobs go wrong, and automatic job re-running if they fail.

To run basespaceapp locally, see the app's help menu!

`basespaceapp --help`


## Other Queue Systems

Basespaceapp is a very modular tool, and is designed to support any execution backend. Right now we only support PBS, but if you know just a little bit of Python, you can add support for any queue easily! *More information coming soon!*

