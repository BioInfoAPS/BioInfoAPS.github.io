---
layout: page
title: Setup
---

{% comment %} Setup {% endcomment %}
<h2 id="setup">Setup</h2>
<p>
  To participate in a this workshop, you will need access to the software described below.
</p>
<p>
  We maintain a list of common issues that occur during installation as a reference for instructors
  that may be useful on the
  <a href = "{{ site.baseSite }}{{ site.faq }}">FAQ page</a>.
</p>

{% comment %} Zoom {% endcomment %}
{% include setup/conference.html %}

{% comment %} SSH Client {% endcomment %}
{% include setup/ssh.html %}

{% comment %} SFTP {% endcomment %}
{% include setup/sftp.html %}

{% comment %} Text editor {% endcomment %}
{% include setup/editor.html %}

{% comment %} Figtree {% endcomment %}
{% include setup/figtree.md %}