from pydantic import BaseModel

from .utils import api_client


class Organization(BaseModel):
    """
    A Rowan organization

    :ivar name: The name of the organization.
    :ivar weekly_credits: The weekly credits of the organization.
    :ivar credits: The credits of the organization.
    """

    name: str
    weekly_credits: float | None = None
    credits: float | None = None


class OrganizationRole(BaseModel):
    """
    A Rowan organization role

    :ivar name: The name of the organization role.
    """

    name: str


class SubscriptionPlan(BaseModel):
    """
    A Rowan subscription plan

    :ivar name: The name of the subscription plan.
    """

    name: str


class IndividualSubscription(BaseModel):
    """
    A Rowan individual subscription

    :ivar subscription_plan: The subscription plan of the individual subscription.
    """

    subscription_plan: SubscriptionPlan


class User(BaseModel):
    """
    A Rowan user

    :ivar uuid: The UUID of the user.
    :ivar username: The username of the user.
    :ivar email: The email of the user.
    :ivar firstname: The first name of the user.
    :ivar lastname: The last name of the user.
    :ivar weekly_credits: The weekly credits of the user.
    :ivar credits: The credits of the user.
    :ivar billing_name: The billing name of the user.
    :ivar billing_address: The billing address of the user.
    :ivar credit_balance_warning: The credit balance warning of the user.
    :ivar organization: The organization of the user.
    :ivar organization_role: The organization role of the user.
    :ivar individual_subscription: The individual subscription of the user.
    """

    uuid: str
    username: str
    email: str
    firstname: str | None = None
    lastname: str | None = None
    weekly_credits: float | None = None
    credits: float | None = None
    billing_name: str | None = None
    billing_address: str | None = None
    credit_balance_warning: float | None = None
    organization: Organization | None = None
    organization_role: OrganizationRole | None = None
    individual_subscription: IndividualSubscription | None = None

    @property
    def name(self) -> str:
        return f"{self.firstname or ''} {self.lastname or ''}".strip()

    @property
    def organization_name(self) -> str | None:
        return self.organization.name if self.organization else None

    @property
    def organization_weekly_credits(self) -> float | None:
        return self.organization.weekly_credits if self.organization else None

    @property
    def organization_credits(self) -> float | None:
        return self.organization.credits if self.organization else None

    @property
    def organization_role_name(self) -> str | None:
        return self.organization_role.name if self.organization_role else None

    @property
    def subscription_plan_name(self) -> str | None:
        return (
            self.individual_subscription.subscription_plan.name
            if self.individual_subscription
            else None
        )

    def credits_available_string(self) -> str:
        """
        Returns a string showing available credits, including organization credits if applicable

        :return: A string showing available credits
        """
        individual_credits = f"Weekly Credits: {self.weekly_credits}\nCredits: {self.credits}"
        if self.organization is not None:
            return (
                individual_credits
                + f"\nOrganization Weekly Credits: {self.organization_weekly_credits}\n"
                f"Organization Credits: {self.organization_credits}"
            )

        return individual_credits


def whoami() -> User:
    """
    Returns the current user
    """

    with api_client() as client:
        response = client.get("/user/me")
        response.raise_for_status()

    return User(**response.json())
