def create_criteria(name, **criteria):
    new_class = type(name, (object,), criteria)
    return new_class