## SAT-and-Reduce Algorithm for Vertex Cover

## Building

Building the program requires Armin Biere's excellent SAT solver [CaDiCaL](https://github.com/arminbiere/cadical) as a dependency. Our Meson-based build system finds CaDiCaL via a pkg-config file. CaDiCaL does not provide such a file, but the following template can be adapted instead:


```
prefix=@PREFIX@

Name: cadical
Version: 1.0.3
Description: CaDiCaL
Cflags: -I\${prefix}/src
Libs: -L\${prefix} -lcadical
```
where `@PREFIX@` is replaced by the build path of CaDiCaL.

With CaDiCaL, `meson` and `ninja` installed, the SAT-and-Reduce algorithm can be build via:

```
mkdir build && cd build
meson --buildtype=debugoptimized ..
ninja
```

and invoked by:

```
./vc-bnb --stats <input file>
```

## Citing

When using this work in an academic context, please cite:

Rick Plachetta and Alexander van der Grinten. SAT-and-Reduce for Vertex Cover: Accelerating Branch-and-Reduce by SAT Solving. ALENEX 2021.

```
@inbook{pg21satreduce,
	author = {Rick Plachetta and Alexander van der Grinten},
	title = {SAT-and-Reduce for Vertex Cover: Accelerating Branch-and-Reduce by SAT Solving},
	booktitle = {2021 Proceedings of the Symposium on Algorithm Engineering and Experiments (ALENEX)},
	chapter = {},
	pages = {169-180},
	doi = {10.1137/1.9781611976472.13}
}
```

## Acknowledgements

We thank Narek Bojikian and Maximilian McKone for work on the implementation of well-known reduction algorithms that our SAT-and-Reduce solver uses as a subprocedure.
