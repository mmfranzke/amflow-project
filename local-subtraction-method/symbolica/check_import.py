from lsmethod.symbolica_backend import check_symbolica_import


def main():
    module = check_symbolica_import()
    print(module)


if __name__ == "__main__":
    main()
