# Folders and projects

Workflows live in a tree of folders, rooted in a project. Folders are how you keep runs grouped and findable. Most of the time `rowan.get_folder("path/to/folder")` is all you need (it creates the path if missing); the rest of this file covers browsing and navigating an existing tree.

## Get a folder by path

```python
folder = rowan.get_folder("CDK2/docking/batch_1")   # creates missing segments by default
folder = rowan.get_folder("CDK2/docking", create=False)  # raises if it doesn't exist
wf = rowan.submit_docking_workflow(..., folder=folder)
```

`get_folder` resolves the path within the active project (see Projects below).

## Browse and traverse an existing tree

Start at the project root, then list or navigate:

```python
root = rowan.root_folder()          # root folder of the active project

root.children()                     # list[Folder]   — subfolders only
root.workflows()                    # list[Workflow] — workflows only
root.contents()                     # list[Folder | Workflow] — both, folders first

docking = root / "CDK2" / "docking" # navigate by name with the / operator
parent = docking.parent()           # Folder, or None at the root
```

- `children()`, `workflows()`, and `contents()` each take a `size=` cap (default 100) and return fully-typed objects you can act on (e.g. `wf.result()`).
- The `/` operator looks up a child folder by exact name. It raises `ValueError` if no child matches, or if **two children share that name** — folder names are not unique within a parent, so in that case use `retrieve_folder(uuid)` with the specific UUID (the error lists the candidates).
- `print_folder_tree()` prints the recursive folder hierarchy (folders only — workflows are not included). Use `children()`/`workflows()`/`contents()` to see workflows at a given level.

## Projects

A project is the top-level container; its root folder is what `root_folder()` returns. The active
project scopes `get_folder` and `root_folder` — set it once, then navigate folders within it.

```python
# Find / list
rowan.list_projects()                       # list[Project]
rowan.list_projects(name_contains="CDK2")   # filtered
proj = rowan.retrieve_project(uuid)         # by UUID (from the project URL)
default = rowan.default_project()           # the account's default project

# Create
proj = rowan.create_project("CDK2 campaign")

# Switch the active project (scopes subsequent get_folder / root_folder)
rowan.set_project("CDK2 campaign")          # by name; sets rowan.project_uuid
rowan.project_uuid = proj.uuid              # or assign the UUID directly

# Rename / delete
proj.update(name="CDK2 campaign v2")
proj.delete()                               # destructive: removes all folders + workflows inside
```

With no project set, the default project is used. `set_project` raises `ValueError` if no project
matches the name exactly.
