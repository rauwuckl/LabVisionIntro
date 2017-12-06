from optparse import OptionParser
import collections
import json


def main():
    parser = OptionParser()
    parser.add_option("-f", "--father", dest="father", help="Father json file. All params not specified in the child file will have this value", action="store")
    parser.add_option("-c", "--child", dest="child", help="Values specified in this file will be used instead of the father ones. (Child values overwrite father values)")
    parser.add_option("-o", "--output", dest="output", help="Name of the output file. If this is not given it will overwrite the child")

    options, args = parser.parse_args()
    is_ok = True
    if options.father is None:
        print("please specify father file")
        is_ok = False
    if options.child is None:
        print("please specify child file")
        is_ok = False

    if not is_ok:
        return

    if options.output is None:
        print("Overwriting the parent file")
        return


    with open(options.father, "r") as f:
        father = json.load(f)

    with open(options.child, "r") as f:
        child = json.load(f)


    new_dict = update(father, child)

    print(new_dict)

    with open(options.output, "w") as f:
        json.dump(new_dict, f, indent=2)




def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


if __name__ == "__main__":
    main()
