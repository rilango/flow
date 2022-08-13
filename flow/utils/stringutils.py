import re


def substitute(input_string, properties):
    """
    Replaces all bind variables with value in properties. If the value is not
    found in properties then it is left without replacement.

    :param input_string: String with binds.
    :param properties: Possible values for the bind variables.
    :return:
    """
    if input_string is None:
        return ''

    matches = re.findall('\{\{[\\w\\-\\.]+\}\}', input_string, re.M | re.I)

    for match in matches:
        bind_variable = match.replace('{{', '').replace('}}', '')
        try:
            value = properties[bind_variable]
            if value is not None:
                input_string = input_string.replace(match, str(value))
        except KeyError as ex:
            pass

    return input_string


def get_bind_variables(input_string):
    """
    Returns all bind variables in the string.

    :param input_string: String with binds.
    :return:
    """
    matches = re.findall('\{\{[\\w\\-\\.]+\}\}', input_string, re.M | re.I)

    bind_variables = []
    for match in matches:
        bind_variable = match.replace('{{', '').replace('}}', '')
        bind_variables.append(bind_variable)
    return bind_variables