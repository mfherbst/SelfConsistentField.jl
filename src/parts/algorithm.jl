"""
    Algorithms need to be setup before usage. Initialization returns a
    ScfIterState which can be iterated.
"""
Base.iterate(::T) where {T <: Algorithm}         = error("Algorithms need to be setup first.")
iterate(::T, ::StepState) where {T <: Algorithm} = error("iterate(::Algorithm, ::Subreport) needs to be implemented by all algorithms.")
