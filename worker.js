// Redirect www.gnomancer.com to the apex host, then serve the static site.
export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    if (url.hostname === "www.gnomancer.com") {
      url.hostname = "gnomancer.com";
      return Response.redirect(url.toString(), 301);
    }
    return env.ASSETS.fetch(request);
  },
};
