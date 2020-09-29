## Gravitational\_Curvatures\_of\_Tesseroids

[![BSD license](http://img.shields.io/badge/license-BSD-lightgrey.svg?style=flat)](https://github.com/xiaoledeng/Gravitational_Curvatures_of_Tesseroids/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/118716751.svg)](https://zenodo.org/badge/latestdoi/118716751)

Developed by [Xiao-Le Deng](http://xiaoledeng.github.io/) in cooperation with [Wen-Bin Shen](http://wbshen.users.sgg.whu.edu.cn/).

Based on the *[Tesseroids](http://tesseroids.leouieda.com/)* software, We offer the codes for modeling **Gravitational Curvatures (GC)** of tesseroids, which are third-order derivatives of the gravitational potential (GP).

>*Tesseroids* is a collection of **command-line tools** for modeling the gravitational potential, acceleration, and gradient (Marussi) tensor. See details at [Leouieda's Github: tesseroids](https://github.com/leouieda/tesseroids).

## Description

`SConstruct` and `src` are the codes for modeling the gravitational curvatures of tesseroids, in corporation with the *[Tesseroids](http://tesseroids.leouieda.com/)* software.

`Experiments` is used to reproduce the results presented in [*Deng and Shen (2019)*](http://dx.doi.org/10.1007/s11200-018-0772-4). 

## License

The codes are used for calculating the gravitational curvatures of a tesseroid under the terms of the
BSD 3-clause license. See [LICENSE](https://github.com/xiaoledeng/Gravitational_Curvatures_of_Tesseroids/blob/master/LICENSE).

## Citing

If you use these codes in your research, please **cite** the papers in your publications:

> Uieda, L., V. Barbosa, and C. Braitenberg (2016), *Tesseroids: Forward-modeling gravitational fields in spherical coordinates*, ***GEOPHYSICS***, F41-F48, doi:[10.1190/geo2015-0204.1](http://dx.doi.org/10.1190/geo2015-0204.1).
> 
> Deng, X.L., Shen, W.B (2019), *Topographic effects up to Gravitational Curvatures of tesseroids: A case study in China*,  ***Studia Geophysica et Geodaetica***, doi:[10.1007/s11200-018-0772-4](http://dx.doi.org/10.1007/s11200-018-0772-4).


## Installing

1. Install the **Tesseroids** software. see help at [leouieda/tesseroids: Forward modeling of gravity fields in spherical coordinates](https://github.com/leouieda/tesseroids);

2. Download the codes at [xiaoledeng/Gravitational\_Curvatures\_of\_Tesseroids](https://github.com/xiaoledeng/Gravitational_Curvatures_of_Tesseroids);

3. Replace the `SConstruct` and put the codes in `src` file at right place, and compile with 

	    $ scons -c
	    $ scons

4. Test the success of installing by using:

		$ tessgzzz -h

		Usage: tessgzzz MODELFILE [OPTIONS]
	
		Calculate the gzzz component due to a tesseroid model on
		specified observation points.
	
		Values are calculated in the local coordinate system of the
		observation point: x-> North  y-> East  z-> Up (away from the
		center of the Earth).
		In order to maintain mainstream convention, component gz is
		calculated with z-> Down.
	
		All units either SI or degrees!
	
		The computation of the gravitational effect of the tesseroids
		is done using the Gauss-Legendre Quadrature (GLQ) numerical
		integration method.
	
		WARNING: Avoid computing directly on top or inside the
	     tesseroids! This will break the GLQ and the formulas!
	
		Input:
  			Computation points passed through standard input (stdin).
  			Reads 3 or more values per line and inteprets the first 3 as:
    		longitude, latitude and height
  			of a computation points. Height should be in meters.
  			Othervalues in the line are ignored.
  			Lines that start with # are ignored as comments.
  			Lines should be no longer than 10000 (ten thousand) characters.

		Output:
  			Printed to standard output (stdout) in the form:
    		lon lat height ... result
  			... represents any values that were read from input and
  			ignored. In other words, the result is appended to the last
  			column of the input. Use this to pipe tessg* programs
  			together.

  			* Comments about the provenance of the data are inserted into

    		the top of the output
    
    	MODELFILE: File containing the tesseroid model
  			* Each tesseroid is specified by the values of its borders

    			and density
  			* The file should contain one tesseroid per line
  			* Each line should have the following column format:

      			West East South North Top Bottom Density
  			* Top and Bottom should be read as 'height to top' and

    			'height to bottom' from the mean Earth radius. Use negative
    			values if bellow the surface, for example when modeling
    			deep structures, and positive if above the surface, for
    			example when modeling topography.
  		* If a line starts with # it will be considered a comment and

    			will be ignored.
    
    	Options:
  			-a             Disable the automatic subdividing of
                 			tesseroids. Subdividing is done to ensure the
                 			GLQ gives accurate results. ONLY USE THIS
                 			OPTION IF YOU KNOW WHAT YOU ARE DOING!
  			-tRATIO        Use a custom distance-size ratio for the
                 			automatic subdivision of tesseroids. ONLY USE
                 			THIS OPTION IF YOU KNOW WHAT YOU ARE DOING!
  			-oOLON/OLAT/OR GLQ order to use in the longitudinal,
                 			latitudinal and radial integrations,
                 			respectively. Defaults to 2/2/2.
                 			Subdividing of tesseroids works best with the
                 			default order.
  			-h             Print instructions.
  			--version      Print version and license information.
  			-v             Enable verbose printing to stderr.
  			-lFILENAME     Print log messages to file FILENAME.

			Part of the Tesseroids package (v2017).
	
			Project site: <http://www.leouieda.com/tesseroids/>
			Report bugs at: <https://github.com/leouieda/tesseroids/issues>
	
			Copyright (C) 2011-2017, Leonardo Uieda.
			This software is distributed under the terms of the BSD License:
			<http://tesseroids.readthedocs.org/en/latest/license.html>
			This is free software: you are free to change and redistribute it.
			There is NO WARRANTY, to the extent permitted by law.
