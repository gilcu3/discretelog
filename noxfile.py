import nox


@nox.session(python=["3.10", "3.11", "3.12", "3.13", "3.14", "pypy3.10", "pypy3.11"])
def tests(session):
    python_executable = session.python
    session.run("poetry", "env", "use", python_executable)
    session.run("poetry", "install", "--with", "dev,test")
    session.run("poetry", "run", "coverage", "run", "-m", "pytest", "--perf")
    session.run("poetry", "run", "coverage", "report")


@nox.session
def lint(session):
    session.run("poetry", "install", "--with", "dev,test")
    session.run("poetry", "run", "black", "--check", ".")
    session.run("poetry", "run", "flake8", ".")


@nox.session
def typing(session):
    session.run("poetry", "install", "--with", "dev,test")
    session.run("poetry", "run", "mypy", ".")
