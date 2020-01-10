from distutils.core import setup, Extension

def main():
    setup(name="RobRepeatFinder",
          version="1.0.0",
          description="Python interface for repeatFinder",
          author="Rob Edwards",
          author_email="raedwards@gmail.com",
          ext_modules=[Extension("RobRepeatFinder", sources=["repeatFinder.cpp"], language='c++')])

if __name__ == "__main__":
    main()
