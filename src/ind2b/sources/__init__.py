"""Standalone clients for the public data sources used by the pipeline.

Each module wraps one provider and returns plain dicts. No client depends on
any hosted agent runtime, so the pipeline runs anywhere Python and network
access are available.
"""
