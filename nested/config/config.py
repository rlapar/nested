import yaml

with open('/etc/nested/config.yml') as stream:
    config = yaml.load(stream)


def args_dict_to_list(args):
    args_list = []
    for arg in args:
        args_list += ['-{}'.format(arg), str(args[arg])]
    return args_list
