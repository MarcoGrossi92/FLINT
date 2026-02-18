# Chemistry Code Generation

FLINT supports both general-purpose chemistry solver and mechanism-specific dedicated routines. The first one is a general subroutine that solves finite-rate chemistry using input data to define all necessary . On the other hand, dedicated routines are source files optimized for a specific chemical mechanism, hard-coding all required info directly in the source files.

The employment of mechanism-specific routines provides significant performance improvements and it is particularly important for production runs.

The generation of a new mechanism routine is always recommended. However, this operation requires the modification of the API, therefore it is suggested for developers or advanced users.

A python program is shipped with FLINT to automatically generate source files starting from yaml files. See the proper [Development](../development/chemistry_generation.md) section for more info.

!!! warning "API modifications required" 
    This procedure generates new FLINT source files and requires rebuilding the library. It is intended for advanced users and developers.

## Performance Comparison

A performance comparison among different methods of integration was performed by testing different chemical mechanisms on 0D batch reactor.

The results clearly show how explicit (hard-coded) routines outperform the other alternatives for each tested mechanism.

Benchmark environment:

- MacBook Air (Apple M1, 8 cores: 4P + 4E)
- 16 GB RAM
- macOS
- Homebrew GCC 15.2.0
- Apple clang 17.0.0 (clang-1700.6.3.2)
- Single-threaded execution

<figure>
  {% include "user/images/barplot.svg" %}
  <figcaption>Normalized execution times</figcaption>
</figure>

