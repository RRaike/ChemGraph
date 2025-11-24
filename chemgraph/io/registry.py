readers = {}
writers = {}


def register_reader(name):
    """Decorator that adds the function to the registry."""

    def decorator(func):
        readers[name] = func
        return func

    return decorator


def register_writer(name):
    """Decorator that adds the function to the registry."""

    def decorator(func):
        writers[name] = func
        return func

    return decorator
