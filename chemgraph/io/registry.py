readers = {}

def register_reader(name):
    """Decorator that adds the function to the registry."""
    def decorator(func):
        readers[name] = func
        return func
    return decorator