Done with the regular course? Here are some more advanced lessons.

{% comment %} List of optional courses {% endcomment %}

<div class="row">
<div class="col-md-6" markdown="1">
<table class="table table-striped"><tbody>

{% for optional in site.optionals %}

<tr><td></td><td><a href="{{ relative_root_path }}{{ optional.url }}">{{ optional.title }}</a></td></tr>

{% endfor %}

</tbody></table>
</div></div>