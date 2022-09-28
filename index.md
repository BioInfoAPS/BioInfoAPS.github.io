---
layout: workshop      # DON'T CHANGE THIS.
# Detailed instructions -> https://carpentries.github.io/workshop-template/customization/index.html
venue: "Plant health 2022"        # brief name of the institution
address: "https://bioinfoaps.github.io/"      # full street address or videoconferencing URL
country: "us"      # lowercase two-letter ISO country code
language: "en"     # lowercase two-letter ISO language code
latitude: "29.637690"        # latitude of workshop venue
longitude: "-82.361710"       # longitude of the workshop venue 
humandate: "Oct 18th, 2022"    # human-readable dates for the workshop 
humantime: "11:00 am - 4:00 pm"    # human-readable times for the workshop 
timezone: "Eastern time"    # human-readable times for the workshop
startdate: 2022-10-18      # machine-readable start date in YYYY-MM-DD
enddate: 2022-10-18        # machine-readable end date in YYYY-MM-DD
instructor: ["Jose C. Huguet-Tapia", "Braham Dhillon"] # boxed, comma-separated list of instructors
moderator: ["Liliana Cano", "Erica M. Goss", "Anuj Sharma"]     # boxed, comma-separated list of helpers
email: ["jhuguet@ufl.edu", "dhillonb@ufl.edu"]    # boxed, comma-separated list of contact email addresses
technical: "anujsharma@ufl.edu"    # email addresss for technical help
collaborative_notes:  # optional: URL for the workshop collaborative notes, e.g. an Etherpad or Google Docs document (e.g., https://pad.carpentries.org/2015-01-01-euphoria)
eventbrite: false     # optional: alphanumeric key for Eventbrite registration, e.g., "1234567890AB" (if Eventbrite is being used)
---

<h2 id="general">General Information</h2>

{% comment %} Introduction {% endcomment %}
{% include intro/intro.md %}

{% comment %} Audience {% endcomment %}
{% include intro/who.md %}

{% comment %} Address {% endcomment %}
{% include intro/where.md %}

{% comment %} Date and Time {% endcomment %}
{: #when}
**When:** {{page.humandate}}, {{page.humantime}} {{page.timezone}}.

{% comment %} Pre-requisite {% endcomment %}
{% include intro/prereq.md %} 

{% comment %} Accessibility {% endcomment %}
{: #accessibility}
**Accessibility:**
We are dedicated to providing a positive and accessible learning environment for all. Please
notify the instructors in advance of the workshop if you require any accommodations or if there is
anything we can do to make this workshop more accessible to you.

{% comment %} Contact {% endcomment %}

[//]: # (Contact)
{: #contact}
**Contact:**
Please email 
{% for email in page.email %}{% if forloop.last and page.email.size > 1 %} or {% else %}{% unless forloop.first %}, {% endunless %}{% endif %} [{{email}}](mailto:{{email}}) {% endfor %}
for more information about the content.
For technical problems regarding the website,
email [{{ page.technical }}](mailto:{{ page.technical }}).

<hr/>


{% comment %} Survey {% endcomment %}
<h2 id="surveys">Surveys</h2>
<p>Please be sure to complete these surveys before and after the workshop.</p>
<p><a href="{{ site.baseSite }}{{ site.survey }}">Pre-workshop Survey</a></p>
<p><a href="{{ site.baseSite }}{{ site.feedback }}">Post-workshop Survey</a></p>
<hr/>

{% comment %} Schedule {% endcomment %}
<h2 id="schedule">Schedule</h2>
{% include intro/schedule.md %}
<hr/>

{% comment %} Optional exercises {% endcomment %}
<h2 id="optional">Optional excercises</h2>
{% include intro/optional_schedule.md %}
<hr/>

{% comment %} Setup {% endcomment %}
<h2 id="setup">Setup</h2>
<p>Go to <a href="setup.html">this page</a> for setup instructions.
