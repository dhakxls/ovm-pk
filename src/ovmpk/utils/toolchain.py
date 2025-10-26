import shutil

def which_all():
    return {
        "smina": shutil.which("smina"),
        "gnina": shutil.which("gnina"),
    }
