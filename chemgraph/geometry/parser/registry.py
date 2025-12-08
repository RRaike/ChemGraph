REGISTRY_GEOMETRY_PARSER = {}


def register_geometry_parser(name):
    """Decorator that adds the function to the registry."""

    def decorator(func):
        if name in REGISTRY_GEOMETRY_PARSER:
            raise ValueError(f"Geometry parser already exists: {name}")

        REGISTRY_GEOMETRY_PARSER[name] = func
        return func

    return decorator
