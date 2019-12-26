from distutils.core import setup, Extension

def main():
    setup(name="repeatfinder",
          version="1.0.0",
          description="Python interface for repeatFinder",
          author="Rob Edwards",
          author_email="raedwards@gmail.com",
          ext_modules=[Extension("repeatFinder", sources=["repeatFinder.cpp"], language='c++')])

if __name__ == "__main__":
    main()
