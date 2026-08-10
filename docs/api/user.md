# User

::: rowan.user
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
      members:
        - User
        - Organization
        - OrganizationRole
        - SubscriptionPlan
        - IndividualSubscription
        - whoami
