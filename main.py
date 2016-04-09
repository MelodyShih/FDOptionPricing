from option import EuropeanCallOption

def main():
	eurocall = EuropeanCallOption()
	eurocall.SetParameter(1.0,0,0,0,0)
	eurocall.SetBC()
	eurocall.SetNumofGrid()
	eurocall.Build()
	eurocall.Solve()

if __name__ == "__main__":
    main()