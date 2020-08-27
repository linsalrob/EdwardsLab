# Django - Getting going

Some brief notes by Rob to start a new web site and application in Django, because you rarely need to start from scratch!

Start a new project in PyCharm, and start a new app in the PyCharm window too.

It should start with `settings.py` and `urls.py` open

In your top level `urls.py` (this should be the file open):
 - add `include` to the django.urls import
 - add a path to your project: `path('passersby/', include('passersby.urls')),` (_Note:_ trailing slash on the path here is _very_ important! [wasted a good hour or two on that!])

In your top level `settings.py`:
 - set time zone to `TIME_ZONE = 'America/Los_Angeles'`
 - under   `INSTALLED_APPS` add your app: `'passersby.apps.PassersbyConfig'`  (If `passersby` is my app name, open `passersby.apps` to find the config method name). PyCharm has likely added that when you added the include above!
 
 
You are now done with those two windows and can close them. Lets focus on the app itself.
 
 
 _ALL OF THIS IS IN THE APP DIRECTORY:_
 
## Create some URLs

Make a new python file called `urls.py` and put this in it. (You can expand this of course):

```
from django.urls import path

from . import views

app_name = 'passersby'
urlpatterns = [
    path('', views.index, name='index'),
]
```

And add an index to `views.py`:
```
from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse


def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")
```





At this point you should be able to open a terminal, and run:

```
python manage.py migrate
```

to update everything and start it up.

You can also test, of course, with 

```
python manage.py runserver
```

 
## Build  a model to get started
 
 See `models.py` for the models set up. 
 
 Once you have created some classes, switch to `appname/admin.py` and
  - add an import e.g. `from .models import Question`
  - register each class with the admin interface: `admin.site.register(Question)`

 
## Set up administrator access
 
 Open a terminal and enter:
 ```
  python manage.py createsuperuser
  ```
  
 
## Set up at least an index view
 
 Remember you can use render to return the object
 Remember to create `templates/appname/` to store the html to render the view
 
# Migrate and test

Now that you have a model and a view, you can run a migrate and test. Don't forget there are two steps, `makemigrations` and then `migrate`. 
 
```
python manage.py makemigrations
python manage.py migrate
python manage.py runserver
```

 
# Install Twitter Bootstrap

```
pip install django-twitter-bootstrap
```

Add `twitter_bootstrap` to `settings.py`:

```
INSTALLED_APPS = [
	...
	'twitter_bootstrap',
]
```




